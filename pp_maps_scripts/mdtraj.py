import os
import subprocess
from glob import glob
from subprocess import PIPE, Popen


def mdtraj_to_pdbs(path_traj: str, path_pdb: str, output_dir: str) -> None:

    if os.path.isdir(output_dir):
        list_files = glob(output_dir + "/*")
        for f in list_files:
            os.remove(f)
    else:
        os.makedirs(output_dir)
    command = [
        "mdconvert",
        path_traj,
        "-t",
        path_pdb,
        "-o",
        output_dir + "/md_.pdb",
    ]
    p = Popen(command, stdin=subprocess.PIPE)
    stdout, stderr = p.communicate()

    with open(output_dir + "/md_.pdb", "r") as all_pdbs_f:
        global_file = all_pdbs_f.readlines()
        list_global = [[]]
        for line in global_file:
            list_global[-1].append(line)
            if line[:6] == "ENDMDL":
                list_global.append([])

    for i in range(len(list_global) - 1):
        with open(f"{output_dir}/md_{i}.pdb", "w") as md_subfile:
            md_subfile.writelines(list_global[i])

    os.remove(output_dir + "/md_.pdb")

