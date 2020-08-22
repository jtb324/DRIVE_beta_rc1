import numpy as np
from numpy.core.numeric import NaN
import pandas as pd
import re
import os.path
from os import path

from graphviz import Digraph
from check_directory import check_dir
from file_exist_checker import Check_File_Exist


class Network_Img_Maker(Check_File_Exist):

    def __init__(self, segments_file, variant_file: str, output_path: str, logger):
        self.file = segments_file
        # This comes from previous singleVariantAnalysis so the file will exist
        self.variant_file_list = variant_file
        self.log_file = logger
        self.output_path = output_path
        self.network_carriers = None
        self.var_of_interest = None

    def isolate_variant_list(self, variant_of_interest: str):
        '''This function will create a list of IIDs who carry a specific variant.'''

        self.var_of_interest = variant_of_interest

        # This uses the reformated variant file
        variant_df = pd.read_csv(self.variant_file_list, sep=",")

        # Isolating the variant_df for the variant of interest
        variant_df_subset = variant_df[variant_df["Variant ID"]
                                       == variant_of_interest]

        iid_list = variant_df_subset["IID"].values.tolist()

        print(f"The number of carriers identified are {len(iid_list)}")

        return iid_list

    def drop_empty_rows(self, loaded_df):
        '''This function just drops empty rows in the dataframe'''
        print(f"dropping empty rows in file {self.file}...")

        nan_value = float("NaN")
        loaded_df.replace("", nan_value, inplace=True)
        loaded_df.dropna(subset=["Pairs"], inplace=True)
        print(loaded_df)
        return loaded_df

    def isolate_ids(self, loaded_df, iid_list):
        '''This function removes the IBD software identifier (The GERMLINE, iLASH, or HapIBD) creates two new rows in the provided dataframe for each set of ids. The function takes the loaded df of individuals as an input'''
        ########################################################
        # Removing the "GERMLINE", "iLASH", and "hapibd" from the Pairs value and then splitting the list:

        self.log_file.info(
            "Creating two new columns in the provide segment file for the two id pairs.")

        loaded_df["Pair_id1"] = loaded_df["Pairs"].str.split(":")

        # Regular expression used to drop GERMLINEq:

        reg_expression = re.compile(
            r'(?=.*GERMLINE)(?=.*iLASH)|(?=.*hapibd)|(?=.*iLASH).*')

        regmatch = np.vectorize(lambda x: bool(reg_expression.match(x)))

        bool_array = regmatch(loaded_df["Pairs"].values)

        loaded_df = loaded_df[bool_array]

        ##########################################################
        # This breaks up the pairs in the dataframe so that they are easily accessible to build networks.
        for index, row in loaded_df.iterrows():

            row['Pair_id1'] = row["Pair_id1"][1]

        loaded_df["Pair_id2"] = loaded_df["Pair_id1"].str.split("-")

        for index, row in loaded_df.iterrows():

            row['Pair_id1'] = row['Pair_id2'][0]
            row['Pair_id2'] = row['Pair_id2'][1]

        # Now have to filter down the file so that all the rows are
        # individuals who are found in the iid_list. Should return a
        # file where Pair_id1 and Pair_id2 are both carriers

        final_df = loaded_df[(loaded_df["Pair_id1"].isin(iid_list)) & (
            loaded_df["Pair_id2"].isin(iid_list))]

        print(
            f"There were {len(final_df)} pairs found within the segment file at {self.file}")
        print("\n")

        self.log_file.info(
            f"There where {final_df} pairs found within the segment file at {self.file}")

        final_df.to_csv(
            "".join([self.output_path, "/", self.var_of_interest, ".pairs.csv"]))

        return final_df
    ################################################################

    def carriers_in_network(self, iid_list, subset_df, ind_in_networks_df):
        '''This function tells the percent of carriers who are in these networks'''

        id1_set = set(subset_df.Pair_id1.values.tolist())

        id2_set = set(subset_df.Pair_id2.values.tolist())

        total_carriers_set = id1_set | id2_set

        carriers_in_network = sum(
            item in total_carriers_set for item in set(iid_list))

        percent_in_networks = carriers_in_network/len(iid_list)*100

        ind_in_networks_df[self.var_of_interest] = percent_in_networks

        # Need to make it so that there is a dataframe with the IID and then a 1 or 0 if it is in a network or not
        # figure out which carriers from the iid_list are in networks
        carriers_in_network_dict = {
            "IID": [],
            "In Network": [],
            "Network ID": []
        }

        # This for loop checks to see if the iid in the iid_list shares a segemnt. It it will then return a 1 if
        # the iid is a carrier or a 0 if it is not avaliable
        for iid in iid_list:

            bool_int = int(iid in total_carriers_set)

            carriers_in_network_dict["IID"].append(iid)
            carriers_in_network_dict["In Network"].append(bool_int)
            carriers_in_network_dict["Network ID"].append(NaN)

        # convert dictionary to Dataframe

        self.network_carriers = pd.DataFrame(
            carriers_in_network_dict, columns=["IID", "In Network", "Network ID"])

        return ind_in_networks_df

    ################################################################

    def draw_networks(self, reformated_df):
        '''This function actually draws the networks. It takes the reformated dataframe from the isolate_ids functions'''

        ##########################################################
        # Logging message
        self.log_file.info(
            "Creating the pdfs of network images and writing them to a directory called network_images")
        ##########################################################
        # Drawing the graph

        # getting a list of all unique IDs in the Pair_id1 column to iterate through
        id_list = reformated_df['Pair_id1'].unique().tolist()

        nodes_visited_set = set()
        edges_drawn = []
        network_number = 1
        # Creating the nodes
        for id1 in id_list:
            print("This is id1")
            print(id1)
            related_graph = Digraph(
                comment="Shared Segment Network", strict=True)

            # Limit the passed dataframe to only those that have a relationship with id1. Only need to Pair_id1 and
            # Pair_id2 columns to make nodes
            nodes_constructed = set()

            if id1 not in nodes_visited_set:

                segments_file_subset = reformated_df[(reformated_df['Pair_id1']
                                                      == id1) | (reformated_df['Pair_id2'] == id1)][["Pair_id1", "Pair_id2"]]

                for row in segments_file_subset.itertuples():

                    Pair_id1 = row[1]
                    Pair_id2 = row[2]
                    print(Pair_id1, Pair_id2)
                    # This if statement only creates the node if it has not been made yet
                    if Pair_id1 not in nodes_visited_set:

                        related_graph.node(Pair_id1, label=Pair_id1)

                        # Keeping track of what nodes have been made
                        nodes_visited_set.add(Pair_id1)

                        nodes_constructed.add(Pair_id1)
                        print(Pair_id1)
                    # This if statement only makes the second node of Pair_id2 if it has not been made yet
                    if Pair_id2 not in nodes_visited_set:

                        related_graph.node(Pair_id2, label=Pair_id2)

                        nodes_visited_set.add(Pair_id2)

                        nodes_constructed.add(Pair_id2)
                        print(Pair_id2)

                    # This checks to see if an edge has been drawn between the two ids
                    if (Pair_id1, Pair_id2) not in edges_drawn:

                        related_graph.edge(Pair_id1, Pair_id2)

                        # Keeping track of edges drawn and the inverse edge
                        edges_drawn.append((Pair_id1, Pair_id2))

                        edges_drawn.append((Pair_id2, Pair_id1))

                    ###########################################################################
                    # Next section is responsible for updating which IID is in which network
                    # Updating the self.network_carriers dataframe so that the Network id gets set to a number

                    idx = self.network_carriers.index[self.network_carriers["IID"] == Pair_id1]

                    self.network_carriers.loc[idx, [
                        "Network ID"]] = str(network_number)

                    idx = self.network_carriers.index[self.network_carriers["IID"] == Pair_id2]

                    self.network_carriers.loc[idx, [
                        "Network ID"]] = str(network_number)

                    ###########################################################################

                # Catching relationships between nodes that are second degree related to id1
                # Need to filter data frame for all the nodes constructed where either Pair_id1 or Pair_id2
                # is in the list but the original id1
                segments_file_subset2 = reformated_df[
                    (reformated_df['Pair_id1'].isin(list(nodes_constructed)) | reformated_df['Pair_id2'].isin(list(nodes_constructed)))]

                segments_file_subset2 = segments_file_subset2[(
                    segments_file_subset2['Pair_id1'] != id1)]

                segments_file_subset2 = segments_file_subset2[(
                    segments_file_subset2['Pair_id2'] != id1)][["Pair_id1", "Pair_id2"]]

                for row in segments_file_subset2.itertuples():

                    Pair_id1 = row[1]
                    Pair_id2 = row[2]

                    # Adding pair_id1 if it has not already been added to both the nodes_visited_set and
                    # nodes_constructed set
                    nodes_visited_set.add(Pair_id1)

                    nodes_constructed.add(Pair_id1)

                    # Then add Pair_id2 to the sets if it is not already in those sets
                    nodes_visited_set.add(Pair_id2)

                    nodes_constructed.add(Pair_id2)

                    print("second ids")
                    print(Pair_id1, Pair_id2)

                    # Checking to see if these second degree edges have
                    # already been drawn. If not the edge is drawn and
                    # then the tuple is appended to the list of edges_drawn

                    if (Pair_id1, Pair_id2) not in edges_drawn:

                        related_graph.edge(Pair_id1, Pair_id2)

                        edges_drawn.append((Pair_id1, Pair_id2))

                        edges_drawn.append((Pair_id2, Pair_id1))

                    # This section accounts for if a total new id was found within the pairs from the segments_file_subset2
                    # This new_id will be Nonetype unless there is a new pair id
                    new_id = None

                    # This sets the new_id value to be Pair_id1 if the node hasn't been visited yet
                    if Pair_id1 not in nodes_constructed:

                        new_id = Pair_id1

                    # If there is a new id then the dataframe has to be subset for the new id
                    # Need to subset the dataframe and then get all the unique values
                    if new_id:
                        print("using new id")
                        # new df subset that contains all variants connected to the new id
                        new_id_segment_subset = reformated_df[
                            (reformated_df['Pair_id1'] == new_id) | (reformated_df["Pair_id2"] == new_id)]

                        # need to get a list of all unique ids in each column
                        uniq_pair_id1 = set(
                            new_id_segment_subset.Pair_id1.values.tolist())

                        uniq_pair_id2 = set(
                            new_id_segment_subset.Pair_id2.values.tolist())

                        uniq_pair_ids = uniq_pair_id1 | uniq_pair_id2

                        uniq_id_subset = reformated_df[(reformated_df['Pair_id1'].isin([uniq_pair_ids]))
                                                       & (reformated_df["Pair_id2"].isin(uniq_pair_ids))]

                        print(uniq_id_subset)

                        for row in uniq_id_subset.itertuples():

                            new_id1 = row[1]
                            new_id2 = row[2]

                            # Adding pair_id1 if it has not already been added to both the nodes_visited_set and
                            # nodes_constructed set
                            nodes_visited_set.add(new_id1)

                            nodes_constructed.add(new_id1)

                            # Then add Pair_id2 to the sets if it is not already in those sets
                            nodes_visited_set.add(new_id2)

                            nodes_constructed.add(new_id2)

                            # Checking to see if these second degree edges have
                            # already been drawn. If not the edge is drawn and
                            # then the tuple is appended to the list of edges_drawn

                            if (new_id1, new_id2) not in edges_drawn:

                                related_graph.edge(new_id1, new_id2)

                                edges_drawn.append((new_id1, new_id2))

                                edges_drawn.append((new_id2, new_id1))

                # This creates a directory for the network pdf files
                img_directory = check_dir(self.output_path, "network_images")

                # making the graphs undirected
                related_graph.edge_attr.update(arrowhead="none")

                # rendering the graphs
                related_graph.render("".join([
                    img_directory, "/", self.var_of_interest, ".network", str(network_number), '.gv']))

                # Clearing the graph for the next loop
                related_graph.clear()

                network_number += 1

                print(nodes_visited_set)

        self.network_carriers.to_csv(
            "".join([self.output_path, "/carriers_in_network", ".", self.var_of_interest, ".csv"]))
