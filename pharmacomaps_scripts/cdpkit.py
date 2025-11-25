import os.path as osp
import sys

import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.MolProp as MolProp
import CDPL.Pharm as Pharm
from rdkit import Chem as RDK_Chem


def split_pdb_receptor_ligand(global_pdb_path: str) -> tuple[str, str]:

    list_ligand = []
    list_receptor = []
    with open(global_pdb_path, "r") as global_pdb:
        list_pdb_lines = global_pdb.readlines()
        for line in list_pdb_lines:
            if line[17:20] == "LIG":
                list_ligand.append(line)
            else:
                list_receptor.append(line)
    complex_path, ext = osp.splitext(global_pdb_path)
    ligand_pdb_path = complex_path + "_ligand" + ext
    receptor_pdb_path = complex_path + "_receptor" + ext
    with open(ligand_pdb_path, "w") as ligand_pdb:
        ligand_pdb.writelines(list_ligand)
    with open(receptor_pdb_path, "w") as receptor_pdb:
        receptor_pdb.writelines(list_receptor)

    return ligand_pdb_path, receptor_pdb_path


def convert_pdb_to_sdf(pdb_path: str) -> str:

    mol = RDK_Chem.MolFromPDBFile(pdb_path, sanitize=True, removeHs=True)
    path, _ = osp.splitext(pdb_path)
    sdf_path = path + ".sdf"
    writer = RDK_Chem.rdmolfiles.SDWriter(sdf_path)
    writer.write(mol)

    return sdf_path


# reads and preprocesses the specified receptor structure
def readAndPrepareReceptorStructure(receptor_file: str) -> Chem.Molecule:
    # create reader for receptor structure (format specified by file extension)
    reader = Chem.MoleculeReader(receptor_file)

    sup_fmts = [
        Chem.DataFormat.MOL2,
        Biomol.DataFormat.PDB,
        Biomol.DataFormat.MMTF,
        Biomol.DataFormat.MMCIF,
    ]

    if (
        reader.getDataFormat() not in sup_fmts
    ):  # check if the format is supported by this script
        sys.exit("Error: receptor input file format not supported")

    rec_mol = (
        Chem.BasicMolecule()
    )  # create an instance of the default implementation of the
    # Chem.Molecule interface that will store the receptor struct.
    try:
        if not reader.read(rec_mol):  # read receptor structure
            sys.exit("Error: reading receptor structure failed")

    except Exception as e:
        sys.exit("Error: reading receptor structure failed:\n" + str(e))

    # preprocess the receptor structure (removal of residues and
    # calculation of properties required by the pharm. generation procedure)
    try:
        rem_atoms = False

        # prepares the receptor structure for pharmacophore generation
        Chem.perceiveSSSR(rec_mol, rem_atoms)
        Chem.setRingFlags(rec_mol, rem_atoms)
        Chem.calcImplicitHydrogenCounts(rec_mol, rem_atoms)
        Chem.perceiveHybridizationStates(rec_mol, rem_atoms)
        Chem.setAromaticityFlags(rec_mol, rem_atoms)

        if Chem.makeHydrogenComplete(
            rec_mol
        ):  # make implicit hydrogens (if any) explicit
            Chem.calcHydrogen3DCoordinates(
                rec_mol
            )  # calculate 3D coordinates for the added expl. hydrogens
            Biomol.setHydrogenResidueSequenceInfo(
                rec_mol, False
            )  # set residue information for the added expl. hydrogens

        MolProp.calcAtomHydrophobicities(
            rec_mol, False
        )  # calculate atom hydrophobicity values (needed for hydrophobic
        # pharm. feature generation)
    except Exception as e:
        sys.exit("Error: processing of receptor structure failed: " + str(e))

    return rec_mol


def pdb_to_pml(receptor_file: str, ligand_file: str, output_pml: str) -> None:

    gen_x_vols = True

    rec_mol = readAndPrepareReceptorStructure(
        receptor_file
    )  # read and preprocess the receptor structure
    lig_reader = Chem.MoleculeReader(
        ligand_file
    )  # create reader for the ligand input file (format specified by file extension)
    ph4_writer = Pharm.FeatureContainerWriter(
        output_pml
    )  # create writer for the generated pharmacophores (format specified by file extension)

    lig_mol = (
        Chem.BasicMolecule()
    )  # create an instance of the default implementation of the
    # Chem.Molecule interface that will store the ligand structures
    ia_ph4 = (
        Pharm.BasicPharmacophore()
    )  # create an instance of the default implementation of the Pharm.Pharmacophore
    # interface that will store the generated pharmacophores

    ph4_gen = (
        Pharm.InteractionPharmacophoreGenerator()
    )  # create an instance of the pharmacophore generator

    ph4_gen.addExclusionVolumes(
        gen_x_vols
    )  # specify whether to generate exclusion volume spheres
    # on pharm. feature atoms of interacting residues
    try:
        i = 1

        # read and process ligand molecules one after the other until the end of input has been reached (or a severe error occurs)
        while lig_reader.read(lig_mol):
            mol_id = Chem.getName(
                lig_mol
            ).strip()  # compose a simple ligand identifier for messages

            if mol_id == "":
                mol_id = "#" + str(
                    i
                )  # fallback if name is empty or not available
            else:
                mol_id = "'%s' (#%s)" % (mol_id, str(i))

            try:
                Pharm.prepareForPharmacophoreGeneration(
                    lig_mol
                )  # make ligand ready for pharm. generation

                ph4_gen.generate(
                    lig_mol, rec_mol, ia_ph4, extract_core_env=True
                )  # generate the pharmacophore (True = extract ligand environment residues on-the-fly)

                try:
                    if not ph4_writer.write(ia_ph4):  # output pharmacophore
                        sys.exit(
                            "Error: writing interaction pharmacophore of molecule %s failed"
                            % mol_id
                        )

                except (
                    Exception
                ) as e:  # handle exception raised in case of severe write errors
                    sys.exit(
                        "Error: writing interaction pharmacophore of molecule %s failed: %s"
                        % (mol_id, str(e))
                    )

            except (
                Exception
            ) as e:  # handle exception raised in case of severe processing errors
                # sys.exit(
                #     "Error: interaction pharmacophore generation for molecule %s failed: %s"
                #     % (mol_id, str(e))
                # )
                raise RuntimeError(
                    "Error: interaction pharmacophore generation for molecule %s failed: %s"
                    % (mol_id, str(e))
                )

            i += 1

    except (
        Exception
    ) as e:  # handle exception raised in case of severe read errors
        raise RuntimeError(
            "Error: reading molecule %s failed: %s" % (str(i), str(e))
        )

    ph4_writer.close()
