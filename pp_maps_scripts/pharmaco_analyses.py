def get_global_dict(list_dict: list[dict]) -> dict:

    dict_global = {}
    for dict_frame in list_dict:
        for interaction_type in dict_frame.keys():
            for interaction_aa in dict_frame[interaction_type].keys():
                interaction = interaction_type + "_" + interaction_aa
                if interaction not in dict_global.keys():
                    dict_global[interaction] = 1
                else:
                    dict_global[interaction] += 1

    return dict_global


def get_interactions_percentage(
    dict_global: dict, number_of_frames: int
) -> dict:

    dict_global_percentage = {}
    for interaction in dict_global.keys():
        dict_global_percentage[interaction] = (
            dict_global[interaction] / number_of_frames * 100
        )

    return dict_global_percentage
