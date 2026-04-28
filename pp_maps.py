import argparse
import configparser
import csv
import json
import multiprocessing as mp
import os
import sys
from glob import glob
from typing import Optional

from biopandas.pdb import PandasPdb

from pp_maps_scripts import cdpkit, ligandscout
from pp_maps_scripts.exceptions import PMLError
from pp_maps_scripts.gmx import traj_to_pdbs
from pp_maps_scripts.map_interactions import generate_heatmap
from pp_maps_scripts.mdtraj import mdtraj_to_pdbs
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
) -> tuple[dict, Optional[dict], Optional[dict]]:
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
    dict_interactions, dict_bb, dict_sc = xml_to_dict(
        pml_path, use_cdpkit, use_ligandscout
    )
    return dict_interactions, dict_bb, dict_sc


def from_traj_to_pharmaco(
    path_traj: str,
    path_pdbs: str,
    path_topol: str,
    output_dir_pdbs: str,
    path_to_pharmacogenerator: str,
    use_ligandscout: bool,
    use_cdpkit: bool,
    use_gromacs: bool,
    use_mdtraj: bool,
    number_processes: int,
) -> tuple[list]:

    if path_traj is not None:
        if use_gromacs:
            traj_to_pdbs(path_traj, path_topol, output_dir_pdbs)

        if use_mdtraj:
            mdtraj_to_pdbs(path_traj, path_topol, output_dir_pdbs)

        list_pdbs = glob(os.path.join(output_dir_pdbs, "*"))
    
    if path_pdbs is not None:
        list_pdbs = glob(os.path.join(path_pdbs, "*"))
    
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
    list_tuple_dict = pool.starmap(loop_core, list_args)

    list_dict = []
    list_bb = []
    list_sc = []
    for tuple_dict in list_tuple_dict:
        list_dict.append(tuple_dict[0])
        list_bb.append(tuple_dict[1])
        list_sc.append(tuple_dict[2])

    return list_dict, list_bb, list_sc


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-xtc",
        type=str,
        help="Path to XTC trajectory file with centered system.",
    )
    parser.add_argument(
        "-pdb",
        type=str,
        help="Path to PDB files directory.",
    )
    parser.add_argument(
        "-topol",
        type=str,
        help="Topology file (TPR for GROMACS or PDB as required for MDTraj).",
        required=True,
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
        "-gromacs",
        action="store_true",
        help="Add this argument to use GROMACS to convert XTC to PDBs.",
    )
    parser.add_argument(
        "-mdtraj",
        action="store_true",
        help="Add this argument to use MDTraj (mdconvert) to convert XTC to PDBs.",
    )
    parser.add_argument(
        "-n",
        type=int,
        default=1,
        help="Number of processors/CPUs to perform the analysis. Default to 1.",
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
    path_pdbs = args.pdb
    path_topol = args.topol
    use_ligandscout = args.ligandscout
    use_cdpkit = args.cdpkit
    use_gromacs = args.gromacs
    use_mdtraj = args.mdtraj
    number_processes = args.n
    if path_traj is not None and path_pdbs is not None:
        print("Too many arguments for input data. Please provide only -xtc or -pdb argument.")
        sys.exit()
    if path_traj is None and path_pdbs is None:
        print("No argument specified for input data. Please provide -xtc or -pdb argument.")
        sys.exit()
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
    if use_gromacs and use_mdtraj:
        print(
            "Too many arguments for trajectory converter. Please select only one tool."
        )
        sys.exit()
    if not use_gromacs and not use_mdtraj:
        print(
            "No argument specified for trajectory converter. Please select one tool."
        )
        sys.exit()
    if use_gromacs and os.path.splitext(path_topol)[1] != ".tpr":
        print("Using GROMACS requires a TPR file as topology. Please provide a .tpr file for argument -topol.")
        sys.exit()
    if use_mdtraj and os.path.splitext(path_topol)[1] != ".pdb":
        print("Using MDTraj requires a PDB file as topology. Please provide a .pdb file for argument -topol.")
        sys.exit()
    output_dir_pdbs = "extracted_pdbs"
    print("Starting PML files creation")
    list_dict_interactions, list_bb, list_sc = from_traj_to_pharmaco(
        path_traj,
        path_pdbs,
        path_topol,
        output_dir_pdbs,
        path_to_pharmacogenerator,
        use_ligandscout,
        use_cdpkit,
        use_gromacs,
        use_mdtraj,
        number_processes,
    )
    global_dict_interactions = get_global_dict(list_dict_interactions)
    global_dict_interactions_percentage = get_interactions_percentage(
        global_dict_interactions, len(list_dict_interactions)
    )
    if use_ligandscout:
        global_dict_bb = get_global_dict(list_bb)
        global_dict_bb_percentage = get_interactions_percentage(
            global_dict_bb, len(list_bb)
        )
        global_dict_sc = get_global_dict(list_sc)
        global_dict_sc_percentage = get_interactions_percentage(
            global_dict_sc, len(list_sc)
        )

    # Generate csv file with interactions and percentages from dict
    with open("pharmacophores.csv", "w") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in global_dict_interactions_percentage.items():
            writer.writerow([key, value])

    if use_ligandscout:
        with open("pharmacophores_bb.csv", "w") as csv_file_bb:
            writer = csv.writer(csv_file_bb)
            for key, value in global_dict_bb_percentage.items():
                writer.writerow([key, value])
        with open("pharmacophores_sc.csv", "w") as csv_file_sc:
            writer = csv.writer(csv_file_sc)
            for key, value in global_dict_sc_percentage.items():
                writer.writerow([key, value])

    # Generate json file with interactions and percentages from dict
    output_json = "pharmacophores.json"
    with open(output_json, "w") as file_interactions:
        json.dump(
            global_dict_interactions_percentage, file_interactions, indent=4
        )

    if use_ligandscout:
        output_json_bb = "pharmacophores_bb.json"
        with open(output_json_bb, "w") as file_interactions_bb:
            json.dump(
                global_dict_bb_percentage, file_interactions_bb, indent=4
            )
        output_json_sc = "pharmacophores_sc.json"
        with open(output_json_sc, "w") as file_interactions_sc:
            json.dump(
                global_dict_sc_percentage, file_interactions_sc, indent=4
            )

    # Heatmap generation from json file
    print("Starting heatmap generation")
    generate_heatmap(output_json, args.output)
    if use_ligandscout:
        generate_heatmap(output_json_bb, "pharmacomap_BB.png")
        generate_heatmap(output_json_sc, "pharmacomap_SC.png")
    print("Pharmacomap successfully generated")
