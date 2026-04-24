import argparse
import collections
import json

import matplotlib.pyplot as plt
import pandas as pd
import seaborn


def sort_amino_acids(aa_list: list) -> tuple[list, dict]:

    dict_aa = {}
    dict_rename = {}

    for aa in aa_list:
        aa_split = aa.split("_")
        dict_aa[int(aa_split[1])] = aa_split[0] + " " + aa_split[1]
        dict_rename[aa] = aa_split[0] + " " + aa_split[1]

    dict_sorted_aa = collections.OrderedDict(sorted(dict_aa.items()))
    return list(dict_sorted_aa.values()), dict_rename


def generate_heatmap(json_file: str, output_file: str) -> None:

    f = open(
        json_file,
    )
    data = json.load(f)
    list_interactions = []
    dict_new_labels = {
        "H": "H",
        "HBD": "D",
        "HBA": "A",
        "PI": "+",
        "NI": "-",
        "AR": "π",
    }

    # Put json file data into a list of interactions list
    for interaction in data.keys():
        list_interactions.append(
            [
                dict_new_labels[interaction.split("_", 1)[0]],
                interaction.split("_", 1)[1],
                data[interaction],
            ]
        )

    # Put list data into a df with 3 colums (interaction, aa, percentage)
    df_interactions = pd.DataFrame(
        list_interactions,
        columns=["Interaction type", "HLA-A2 amino acid", "Frames percentage"],
    )

    # Set representation of data (x and y axes, heatmap data)
    tableau = df_interactions.pivot(
        index="Interaction type",
        columns="HLA-A2 amino acid",
        values="Frames percentage",
    )
    ordered_aa_list, dict_rename = sort_amino_acids(tableau.columns.to_list())
    tableau = tableau.rename(columns=dict_rename)
    tableau = tableau[ordered_aa_list]

    # Set figure size
    fig, ax = plt.subplots(figsize=(20, 4))

    # Generate heatmap from df
    map_interactions = seaborn.heatmap(
        tableau,
        annot=True,
        fmt=".0f",
        ax=ax,
        cmap=seaborn.cubehelix_palette(as_cmap=True),
    )

    # Add title to the figure according to json file name
    # ax.set_title(
    #    "Pharmacophore analysis (with percentage of frames) of "
    #    + json_file.split(".")[0]
    #    + " complex"
    # )
    interactions_fig = map_interactions.get_figure()

    # Save heatmap into file with adjustment of image size
    interactions_fig.savefig(output_file, bbox_inches="tight")
