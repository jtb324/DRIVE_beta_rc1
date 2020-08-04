import numpy as np
from graphviz import Digraph
from check_directory import check_dir
from file_exist_checker import Check_File_Exist


class Network_Img_Maker(Check_File_Exist):

    def __init__(self, segments_file, output_path, logger):
        self.file = segments_file
        self.log_file = logger
        self.output_path = output_path

    def isolate_ids(self, loaded_df):
        '''This function removes the IBD software identifier (The GERMLINE, iLASH, or HapIBD) creates two new rows in the provided dataframe for each set of ids. The function takes the loaded df of individuals as an input'''
        ########################################################
        # Removing the "GERMLINE", "iLASH", and "hapibd" from the Pairs value and then splitting the list:

        self.log_file.info(
            "Creating two new columns in the provide segment file for the two id pairs.")

        loaded_df["Pair_id1"] = loaded_df["Pairs"].str.split(":")

        for index, row in loaded_df.iterrows():

            row['Pair_id1'] = row["Pair_id1"][1]

        loaded_df["Pair_id2"] = loaded_df["Pair_id1"].str.split("-")

        for index, row in loaded_df.iterrows():

            row['Pair_id1'] = row['Pair_id2'][0]
            row['Pair_id2'] = row['Pair_id2'][1]

        return loaded_df

    def draw_networks(self, reformated_df):
        '''This function actually draws the networks. It takes the reformated dataframe from the isolate_ids functions'''

        ##########################################################
        # Drawing the graph
        related_graph = Digraph(comment="Shared Segment Network", strict=True)

        # getting a list of all unique IDs in the Pair_id1 column to iterate through
        id_list = reformated_df['Pair_id1'].unique().tolist()

        nodes_visited = []

        # Creating the nodes
        for id1 in id_list:
            current_node_name = id1

            if id1 not in nodes_visited:

                related_graph.node(id1, label=id1)

            segments_file_subset = reformated_df[(reformated_df['Pair_id1']
                                                  == id1) | (reformated_df['Pair_id2'] == id1)]

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

            segments_file_subset2 = reformated_df[(
                reformated_df['Pair_id1'].isin(unique_id_list)) & (reformated_df['Pair_id2'].isin(id2_list))]

            for row in segments_file_subset2[['Pair_id1', 'Pair_id2']].itertuples():

                if row[1] not in nodes_visited or row[2] not in nodes_visited:

                    related_graph.node(row[1], label=row[1])

                    related_graph.node(row[2], label=row[2])

                related_graph.edge(row[1], row[2])

            img_directory = check_dir(self.output_path, "network_images")

            related_graph.render(img_directory+"/"+id1+'.gv')
