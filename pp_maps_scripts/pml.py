import xmltodict

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


def parse_xml(data_path: str) -> dict:

    with open(data_path) as xml_file:
        xml_data = xml_file.read()

    # Convert XML to Python dictionary
    dict_data = xmltodict.parse(xml_data)

    return dict_data


def xml_to_dict(data_path: str, use_cdpkit: bool) -> dict[str, dict]:

    dict_data = parse_xml(data_path)
    if use_cdpkit:
        dict_data = dict_data["ElementContainer"]["ContainerPharmacophores"][
            "alignmentElement"
        ]
    dict_res_p = count_interactions("point", dict_data)
    dict_res_v = count_interactions("vector", dict_data)

    dict_res_p_v = {**dict_res_p, **dict_res_v}

    return dict_res_p_v
