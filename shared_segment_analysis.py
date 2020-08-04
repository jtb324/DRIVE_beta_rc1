# %%
import pandas as pd
import numpy as np
from graphviz import Digraph
from check_directory import check_dir

########################################################
# read in the shared segment file for a specific variant
# %%
segments_file = pd.read_csv(
    "../shared_segments/reformat_IBD_exm652227.116988685_117149147.allpair.txt", sep=",", skiprows=5, header=None)

########################################################
# Dropping unnecessary index file and then renaming the column
segments_file = segments_file.drop([0], axis=1).rename(columns={1: "Pairs"})

########################################################
# Removing the "GERMLINE", "iLASH", and "hapibd" from the Pairs value and then splitting the list:
segments_file["Pair_id1"] = segments_file["Pairs"].str.split(":")
# %%
for index, row in segments_file.iterrows():

    row['Pair_id1'] = row["Pair_id1"][1]

segments_file["Pair_id2"] = segments_file["Pair_id1"].str.split("-")

for index, row in segments_file.iterrows():

    row['Pair_id1'] = row['Pair_id2'][0]
    row['Pair_id2'] = row['Pair_id2'][1]
# %%
##########################################################
# Drawing the graph
related_graph = Digraph(comment="Shared Segment Network", strict=True)

# getting a list of all unique IDs in the Pair_id1 column to iterate through
id_list = segments_file['Pair_id1'].unique().tolist()

nodes_visited = []

# Creating the nodes
for id1 in id_list:
    current_node_name = id1

    if id1 not in nodes_visited:

        related_graph.node(id1, label=id1)

    segments_file_subset = segments_file[(segments_file['Pair_id1']
                                          == id1) | (segments_file['Pair_id2'] == id1)]

    id2_list = segments_file_subset["Pair_id2"].unique().tolist()

    unique_id_list = np.setdiff1d(
        segments_file_subset["Pair_id1"], id2_list).tolist()

    nodes_visited.append(id1)

    for id2 in id2_list:

        if id2 not in nodes_visited and id2 != id1:

            related_graph.node(id2, label=id2)

        related_graph.edge(id1, id2)

        nodes_visited.append(id2)

    for variant_id in unique_id_list:

        related_graph.node(variant_id, label=variant_id)

        related_graph.edge(variant_id, id1)

        nodes_visited.append(variant_id)

    # This section will try to match up the ids in unique_id_list where
    # the id in Pair_id2 column does not match id1

    segments_file_subset2 = segments_file[(
        segments_file['Pair_id1'].isin(unique_id_list)) & (segments_file['Pair_id2'].isin(id2_list))]

    for row in segments_file_subset2[['Pair_id1', 'Pair_id2']].itertuples():
        if row[1] not in nodes_visited or row[2] not in nodes_visited:

            related_graph.node(row[1], label=row[1])

            related_graph.node(row[2], label=row[2])

        related_graph.edge(row[1], row[2])

    img_directory = check_dir("./", "network_images")

    related_graph.render(img_directory+"/"+id1+'.gv')


# %%
