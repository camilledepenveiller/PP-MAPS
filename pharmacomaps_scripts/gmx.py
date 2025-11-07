import os
import subprocess
from glob import glob
from subprocess import PIPE, Popen
from time import sleep


def traj_to_pdbs(path_traj: str, path_tpr: str, output_dir: str) -> None:

    if os.path.isdir(output_dir):
        list_files = glob(output_dir + "/*")
        for f in list_files:
            os.remove(f)
    else:
        os.makedirs(output_dir)
    command = [
        "gmx",
        "trjconv",
        "-f",
        path_traj,
        "-s",
        path_tpr,
        "-sep",
        "-o",
        output_dir + "/md_.pdb",
    ]
    p = Popen(command, stdin=subprocess.PIPE)
    sleep(1)
    p.communicate(input=b"1")
    stdout, stderr = p.communicate()
