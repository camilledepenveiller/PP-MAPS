from typing import Optional

import xmltodict
from biopandas.pdb import PandasPdb

# To do if needed : pip install xmltodict


def count_interactions(category: str, dict_data: dict) -> dict:

    dict_res = {}

    for i in range(len(dict_data["pharmacophore"][category])):
        interaction_type = dict_data["pharmacophore"][category][i]["@name"]
        if interaction_type not in dict_res.keys():
            dict_res[interaction_type] = {}
        aa_list = dict_data["pharmacophore"][category][i][
            "@envCompound"
        ].split(",")
        for a in range(len(aa_list)):
            aa_list[a] = aa_list[a].strip()
            if aa_list[a] not in dict_res[interaction_type].keys():
                dict_res[interaction_type][aa_list[a]] = 1
            else:
                dict_res[interaction_type][aa_list[a]] += 1

    return dict_res


def count_interactions_BB_SC(
    category: str, dict_data: dict
) -> tuple[dict, dict]:

    ppdb = PandasPdb().read_pdb("extracted_pdbs/md_0.pdb")
    pdb_df = ppdb.df["ATOM"]
    list_BB_atoms = [
        "N",
        "H",
        "H1",
        "H2",
        "H3",
        "CA",
        "HA1",
        "HA",
        "C",
        "O",
        "OC1",
        "OC2",
    ]
    dict_bb = {}
    dict_sc = {}

    # Browse PML file and extracts data for each interaction type from category
    for i in range(len(dict_data["pharmacophore"][category])):
        interaction_type = dict_data["pharmacophore"][category][i]["@name"]
        if interaction_type not in dict_bb.keys():
            dict_bb[interaction_type] = {}
        if interaction_type not in dict_sc.keys():
            dict_sc[interaction_type] = {}

        # Create a list of receptor atoms (j) implied in interactions
        list_atoms_env = []
        if isinstance(
            dict_data["pharmacophore"][category][i]["environmentAtom"], list
        ):
            for j in range(
                len(dict_data["pharmacophore"][category][i]["environmentAtom"])
            ):
                list_atoms_env.append(
                    dict_data["pharmacophore"][category][i]["environmentAtom"][
                        j
                    ]["@reference"]
                )
        else:
            list_atoms_env.append(
                dict_data["pharmacophore"][category][i]["environmentAtom"][
                    "@reference"
                ]
            )

        # Create a list of interactions of one type with residue id and aa atom type (BB or SC)
        list_interactions_occurrences = []
        for atom in list_atoms_env:
            atom_type = pdb_df["atom_name"][int(atom) - 1]
            aa_name = pdb_df["residue_name"][int(atom) - 1]
            aa_number = pdb_df["residue_number"][int(atom) - 1]
            aa_chain = pdb_df["chain_id"][int(atom) - 1]
            if atom_type in list_BB_atoms:
                aa_bb_sc = "BB"
            else:
                aa_bb_sc = "SC"
            interaction_id = (
                aa_name
                + "_"
                + str(aa_number)
                + "_"
                + aa_chain
                + "_"
                + aa_bb_sc
            )
            if interaction_id not in list_interactions_occurrences:
                list_interactions_occurrences.append(interaction_id)

        # Fill dictionaries bb and sc with corresponding interactions occurrences for one type
        for interaction in list_interactions_occurrences:
            if interaction.endswith("BB"):
                if interaction not in dict_bb[interaction_type].keys():
                    dict_bb[interaction_type][interaction] = 1
                else:
                    dict_bb[interaction_type][interaction] += 1
            else:
                if interaction not in dict_sc[interaction_type].keys():
                    dict_sc[interaction_type][interaction] = 1
                else:
                    dict_sc[interaction_type][interaction] += 1

    return dict_bb, dict_sc


def parse_xml(data_path: str) -> dict:

    with open(data_path) as xml_file:
        xml_data = xml_file.read()

    # Convert XML to Python dictionary
    dict_data = xmltodict.parse(xml_data)

    return dict_data


def xml_to_dict(
    data_path: str, use_cdpkit: bool, use_ligandscout: bool
) -> tuple[dict, Optional[dict], Optional[dict]]:

    dict_data = parse_xml(data_path)
    if use_cdpkit:
        dict_data = dict_data["ElementContainer"]["ContainerPharmacophores"][
            "alignmentElement"
        ]
    dict_res_p = count_interactions("point", dict_data)
    dict_res_v = count_interactions("vector", dict_data)

    dict_res_p_v = {**dict_res_p, **dict_res_v}
    dict_res_p_v_bb = None
    dict_res_p_v_sc = None

    if use_ligandscout:
        dict_res_p_bb, dict_res_p_sc = count_interactions_BB_SC(
            "point", dict_data
        )
        dict_res_v_bb, dict_res_v_sc = count_interactions_BB_SC(
            "vector", dict_data
        )

        dict_res_p_v_bb = {**dict_res_p_bb, **dict_res_v_bb}
        dict_res_p_v_sc = {**dict_res_p_sc, **dict_res_v_sc}

    return dict_res_p_v, dict_res_p_v_bb, dict_res_p_v_sc


# if __name__ == "__main__":

#     dict_test = parse_xml("tmp/output_SpawnPoolWorker-1.pml")
#     count_interactions_BB_SC("point", dict_test)
