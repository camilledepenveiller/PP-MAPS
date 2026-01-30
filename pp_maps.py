import argparse
import configparser
import json
import multiprocessing as mp
import os
import sys
from glob import glob

from biopandas.pdb import PandasPdb

from pp_maps_scripts import cdpkit, ligandscout
from pp_maps_scripts.exceptions import PMLError
from pp_maps_scripts.gmx import traj_to_pdbs
from pp_maps_scripts.map_interactions import generate_heatmap
from pp_maps_scripts.pdb_modif import modify_pdb
from pp_maps_scripts.pharmaco_analyses import (
    get_global_dict,
    get_interactions_percentage,
)
from pp_maps_scripts.pml import xml_to_dict

PML_TMP_DIRECTORY = "tmp"


def loop_core(
    pdb_path: str,
    path_to_pharmacogenerator: str,
    use_ligandscout: bool,
    use_cdpkit: bool,
) -> dict[str, dict]:
    ppdb = PandasPdb().read_pdb(pdb_path)
    ppdb_modified = modify_pdb(ppdb)
    ppdb_modified.to_pdb(pdb_path)

    pml_name = "output_" + mp.current_process().name + ".pml"
    pml_path = os.path.join(PML_TMP_DIRECTORY, pml_name)

    if use_ligandscout:
        ligandscout.pdb_to_pml(pdb_path, path_to_pharmacogenerator, pml_path)

    if use_cdpkit:
        ligand_path, receptor_path = cdpkit.split_pdb_receptor_ligand(pdb_path)
        sdf_path = cdpkit.convert_pdb_to_sdf(ligand_path)
        cdpkit.pdb_to_pml(receptor_path, sdf_path, pml_path)

    if not os.path.isfile(pml_path):
        raise PMLError(
            message=f"PML path ({pml_path}) not found. "
            "Check your pharmacophore generator installation/license."
        )
    dict_interactions = xml_to_dict(pml_path, use_cdpkit)
    return dict_interactions


def from_traj_to_pharmaco(
    path_traj: str,
    path_tpr: str,
    output_dir_pdbs: str,
    path_to_pharmacogenerator: str,
    use_ligandscout: bool,
    use_cdpkit: bool,
    number_processes: int,
) -> list[dict]:

    traj_to_pdbs(path_traj, path_tpr, output_dir_pdbs)

    list_pdbs = glob(output_dir_pdbs + "/*")
    list_args = [
        (pdb, path_to_pharmacogenerator, use_ligandscout, use_cdpkit)
        for pdb in list_pdbs
    ]

    if os.path.isdir(PML_TMP_DIRECTORY):
        list_files = glob(os.path.join(PML_TMP_DIRECTORY, "*"))
        for f in list_files:
            os.remove(f)
    else:
        os.makedirs(PML_TMP_DIRECTORY)

    pool = mp.Pool(number_processes)
    list_dict = pool.starmap(loop_core, list_args)

    return list_dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-xtc",
        type=str,
        help="XTC trajectory file with centered system.",
        required=True,
    )
    parser.add_argument(
        "-tpr", type=str, help="TPR file used for MD.", required=True
    )
    parser.add_argument(
        "-ligandscout",
        action="store_true",
        help="Add this argument to use LigandScout as pharmacophore generator.",
    )
    parser.add_argument(
        "-cdpkit",
        action="store_true",
        help="Add this argument to use CDPKit as pharmacophore generator.",
    )
    parser.add_argument(
        "-n",
        type=int,
        default=1,
        help="Number of processes to perform the analysis. Default to 1.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="PNG output file with heatmap of pharmacophore features of the whole MD. Default to pharmacomap.png",
        default="pharmacomap.png",
    )
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read("config.ini")
    path_to_pharmacogenerator = config["PHARMACO_CONFIG"][
        "Pharmacophore_generator_path"
    ]

    path_traj = args.xtc
    path_tpr = args.tpr
    use_ligandscout = args.ligandscout
    use_cdpkit = args.cdpkit
    number_processes = args.n
    if use_ligandscout and use_cdpkit:
        print(
            "Too many arguments for pharmacophore generator. Please select only one tool."
        )
        sys.exit()
    if not use_ligandscout and not use_cdpkit:
        print(
            "No argument specified for pharmacophore generator. Please select one tool."
        )
        sys.exit()
    output_dir_pdbs = "extracted_pdbs"
    print("Starting PML files creation")
    list_dict_interactions = from_traj_to_pharmaco(
        path_traj,
        path_tpr,
        output_dir_pdbs,
        path_to_pharmacogenerator,
        use_ligandscout,
        use_cdpkit,
        number_processes,
    )
    global_dict_interactions = get_global_dict(list_dict_interactions)
    global_dict_interactions_percentage = get_interactions_percentage(
        global_dict_interactions, len(list_dict_interactions)
    )
    output_json = "pharmacophores.json"
    with open(output_json, "w") as file_interactions:
        json.dump(
            global_dict_interactions_percentage, file_interactions, indent=4
        )
    print("Starting heatmap generation")
    generate_heatmap(output_json, args.output)
    print("Pharmacomap successfully generated")
