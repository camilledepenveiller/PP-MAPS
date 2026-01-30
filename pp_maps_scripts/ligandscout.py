import os.path as osp
from subprocess import PIPE, Popen


def pdb_to_pml(
    pdb_file: str, path_to_ipharmgen: str, output_path: str
) -> None:

    ipharmgen_path = osp.join(path_to_ipharmgen, "ipharmgen")
    command = [
        ipharmgen_path,
        "-i",
        pdb_file,
        "-o",
        output_path,
    ]
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
