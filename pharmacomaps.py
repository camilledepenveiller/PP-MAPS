import argparse
import configparser
import json
import multiprocessing as mp
import os
from glob import glob

from biopandas.pdb import PandasPdb

from pharmacomaps_scripts.gmx import traj_to_pdbs
from pharmacomaps_scripts.ligandscout import pdb_to_pml
from pharmacomaps_scripts.map_interactions import generate_heatmap
from pharmacomaps_scripts.pdb_modif import modify_pdb
from pharmacomaps_scripts.pharmaco_analyses import (
    get_global_dict,
    get_interactions_percentage,
)
from pharmacomaps_scripts.pml import xml_to_dict

PML_TMP_DIRECTORY = "tmp"


def loop_core(pdb_path: str, path_to_ipharmgen: str) -> dict[str, dict]:
    ppdb = PandasPdb().read_pdb(pdb_path)
    ppdb_modified = modify_pdb(ppdb)
    ppdb_modified.to_pdb(pdb_path)

    pml_name = "output_" + str(mp.current_process().name) + ".pml"
    pml_path = os.path.join(PML_TMP_DIRECTORY, pml_name)

    pdb_to_pml(pdb_path, path_to_ipharmgen, pml_path)

    dict_interactions = xml_to_dict(pml_path)
    return dict_interactions


def from_traj_to_pharmaco(
    path_traj: str,
    path_tpr: str,
    output_dir_pdbs: str,
    path_to_ipharmgen: str,
    number_processes: int,
) -> list[dict]:

    traj_to_pdbs(path_traj, path_tpr, output_dir_pdbs)

    list_pdbs = glob(output_dir_pdbs + "/*")
    list_args = [(pdb, path_to_ipharmgen) for pdb in list_pdbs]

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
    path_to_ipharmgen = config["PHARMACO_CONFIG"][
        "Pharmacophore_generator_path"
    ]

    path_traj = args.xtc
    path_tpr = args.tpr
    number_processes = args.n
    output_dir_pdbs = "extracted_pdbs"
    print("Starting PML files creation")
    list_dict_interactions = from_traj_to_pharmaco(
        path_traj,
        path_tpr,
        output_dir_pdbs,
        path_to_ipharmgen,
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
