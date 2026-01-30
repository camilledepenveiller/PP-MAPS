import os.path as osp
from subprocess import PIPE, Popen

from pp_maps_scripts.exceptions import IPHARMGENError


def pdb_to_pml(
    pdb_file: str, path_to_ipharmgen: str, output_path: str
) -> None:

    ipharmgen_path = osp.join(path_to_ipharmgen, "ipharmgen")
    if not osp.isfile(ipharmgen_path):
        raise IPHARMGENError(
            message=f"IPHARMGEN path ({ipharmgen_path}) not found. "
        )
    command = [
        ipharmgen_path,
        "-i",
        pdb_file,
        "-o",
        output_path,
    ]
    p = Popen(command, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
