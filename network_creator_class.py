import numpy as np
from numpy.core.numeric import NaN
import pandas as pd
import re
import os.path
from os import path
from graphviz import Digraph
import glob

from check_directory import check_dir
from file_exist_checker import Check_File_Exist


class Network_Img_Maker(Check_File_Exist):

    def __init__(self, segments_file, variant_file: str, output_path: str, logger):
        # This will be the directory to where all the allpair.new.txt files are
        self.file = segments_file
        # This comes from previous singleVariantAnalysis so the file will exist
        self.variant_file_list = variant_file
        self.log_file = logger
        self.output_path = output_path
        self.network_carriers = None
        self.var_of_interest = None
        # This is used in creating an id set for drawing the networks later
        self.id_list = None
        self.curr_dir = os.getcwd()

    def gather_allpair_files(self, file_directory: str, file_tag: str):
        '''This function will gather all of the allpair.new.txt files which contain information about pairs. It will also be used to get the 'chr#_list.single_variant.csv' files.'''
        os.chdir(file_directory)

        file_list = []

        for file in glob.glob(file_tag):

            full_file_path = "".join([file_directory, file])

            print(full_file_path)

            file_list.append(full_file_path)

        os.chdir(self.curr_dir)

        return file_list

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

                ids1_set = set(segments_file_subset.Pair_id1.values.tolist())

                ids2_set = set(segments_file_subset.Pair_id2.values.tolist())

                ids_list = ids1_set | ids2_set

                # insert recursive function here
                self.generate_pair_list(
                    ids_list, reformated_df, id1)

                print(self.id_list)

                # Subsetting the original dataframe for this full list
                full_subset_df = reformated_df[(reformated_df.Pair_id1.isin(self.id_list)) | (
                    reformated_df.Pair_id2.isin(self.id_list))][["Pair_id1", "Pair_id2"]]

                # iterating through each row
                for row in full_subset_df.itertuples():

                    # breaking up the tuple into each pair
                    Pair_id1 = row[1]
                    Pair_id2 = row[2]

                    # This if statements only makes Pair_id1 a node if it is not already constructed
                    if Pair_id1 not in nodes_constructed:

                        related_graph.node(Pair_id1, label=Pair_id1)

                        # Keeping track of what nodes have been made
                        nodes_visited_set.add(Pair_id1)

                        nodes_constructed.add(Pair_id1)

                    # This if statement only makes the second node of Pair_id2 if it has not been made yet
                    if Pair_id2 not in nodes_visited_set:

                        related_graph.node(Pair_id2, label=Pair_id2)

                        # Keeps track of what nodes have been made
                        nodes_visited_set.add(Pair_id2)

                        nodes_constructed.add(Pair_id2)

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

                # This creates a directory for the network pdf files
                img_directory = check_dir(self.output_path, "network_images")

                # making the graphs undirected
                related_graph.edge_attr.update(arrowhead="none")

                # rendering the graphs
                related_graph.render("".join([
                    img_directory, "/", self.var_of_interest, ".network", str(network_number), '.gv']))

                # Clearing the graph for the next loop
                related_graph.clear()

                # This updates which network is being created
                network_number += 1

        # This writes all the network_carriers df to a csv file
        self.network_carriers.to_csv(
            "".join([self.output_path, "/carriers_in_network", ".", self.var_of_interest, ".csv"]))

    def generate_pair_list(self, current_id_set: set, segment_df: pd.DataFrame, current_id1: str) -> set:
        '''This function identifies all the potential pair ids for '''

        # First subset the dataframe for all the current ids
        new_df_subset = segment_df[(segment_df['Pair_id1'].isin(current_id_set)) | (
            segment_df['Pair_id2'].isin(current_id_set))][["Pair_id1", "Pair_id2"]]
        print("Current id set")
        print(current_id_set)
        # These next three lines use sets to avoid duplicate values
        # This gets all ids from the pair one column
        id1_set = set(new_df_subset.Pair_id1.values.tolist())
        # This gets all ids from the Pair two column
        id2_set = set(new_df_subset.Pair_id2.values.tolist())
        # this line takes teh union of each set to get all unique values
        total_id_set = id1_set | id2_set
        print("total id set")
        print(total_id_set)

        # If the total_id_set equals the current_id_set then the function stops
        if total_id_set == current_id_set:
            print(f"found all ids connected to the probe {current_id1}")
            # returns the set of ids
            self.id_list = total_id_set
            return
        # If the two sets are not equal then it recursively calls the function again
        else:

            self.generate_pair_list(total_id_set, segment_df, current_id1)
