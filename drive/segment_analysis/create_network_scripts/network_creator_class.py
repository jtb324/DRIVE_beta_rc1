import numpy as np
from numpy.core.numeric import NaN
import pandas as pd
import os.path
from os import path
from typing import List, Dict

import utility_scripts
class Network_Prep:
    
    def __init__(self, allpair_file_dir: str, network_dir: str, analysis_type: str) -> None:
        self.allpair_file_list: list = utility_scripts.get_file_list(allpair_file_dir, "*allpair.txt")
        self.output: str = network_dir
        self.analysis_type: str = analysis_type


    def determine_carriers(self, confirmed_carriers_file: str):
        """Function to get a dictionary of all of the carriers per variant/gene per chromosome
        Parameters
        __________
        confirmed_carriers_file : str
            string that list the filepath to the confirmed_carriers.txt file that has the iids and if they are a confirmed carrier of a variant or a gene.
        
        Returns 
        _______
        Dict[str, Dict[str, List[str]]]
            returns a dictionary where the key is the chromosome and the values is an inner dictionary where the key is the variant or gene and the value is a list of iids that carry that variant
        """
        # creating a dictionary that can keep list of suspected carriers
        iid_dict: Dict[str, Dict[str, List[str]]] = {}
        
        # loading the carriers file into a dataframe
        carriers_df: pd.DataFrame = pd.read_csv(confirmed_carriers_file, sep="\t")

        # if the analysis type is phenotype then it will make the inner key the 
        # gene name else it will make the inner key the variant id
        if self.analysis_type == "phenotype":
            
            # getting a list of all the genes from the carriers_df

            gene_list: List[str] = list(set(carriers_df["gene_name"].values.tolist()))
        
            for gene in gene_list:

                # getting the chromosome number for the specific variant
                chr_num: str = str(list(set(carriers_df[carriers_df["gene_name"] == gene].chr.values.tolist()))[0])

                if len(chr_num) == 1:
                    chr_num = "".join(["0", chr_num])

                # getting the list of iids associated with that variant and ch    romosome
                iid_list: List[str] = carriers_df[carriers_df["gene_name"] == gene].IID.values.tolist()

                # writing this list ot a dictionary
                iid_dict["".join(["chr",chr_num])] = {gene: iid_list}

        else:   
            # getting a list of all the variants from the carriers_df

            variant_list: List[str] = list(set(carriers_df["variant_id"].values.tolist()))

            for variant in variant_list:

                # getting the chromosome number for the specific variant
                chr_num: str = list(set(carriers_df[carriers_df["variant_id"] == variant].chr.values.tolist()))[0]

                # getting the list of iids associated with that variant and chromosome
                iid_list: List[str] = carriers_df[carriers_df["variant_id"] == variant].IID.values.tolist()

                # set a default value for the dictionary 
                iid_dict.setdefault("".join(["chr",str(chr_num)]), {})

                # writing this list ot a dictionary
                iid_dict["".join(["chr",str(chr_num)])][variant] = iid_list

        # creating an attribute that keeps track of these values
        self.iid_dict : Dict[str, Dict[str, List[str]]] = iid_dict


# need a class to keep track of the networks
class Network_Maker:
    """class to determine the different networks"""

    def __init__(self, chr_num: str, identifier: str, iid_list: List[str], analysis_type: str) -> None:
        self.chr_num: str = chr_num
        self.identifier: str = identifier
        self.iid_list: List[str] = iid_list
        self.analysis_type: str = analysis_type
    
    @staticmethod
    def fix_chr(chr_number: str) -> str:
        """Function to make sure that the chromosome number if it is of the form chrX gets converted to 
        chrXX and returns that and if the format is already correct then it just returns the string. The
        X refers to a digit
        Parameter
        _________
        chr_number : str
            string that has the chromosome number
        
        Returns 
        _______
        str
            returns the properly formatted chromosome number string
        """
        if len(chr_number) == 4:
            return "".join([chr_number[:3], "0", chr_number[-1]])

        else:
            return chr_number

    def find_allpair_file(self, allpair_list: List[str]) -> str:
        """Function to find the specific allpair file that lines up with the chr number and the variant/gene name
        Parameters
        __________
        allpair_list : List[str]
            list of all the allpair files
        
        Returns
        _______
        str
            specific filepath to the allpair file that matches the 
            chr_num and variant_id
        """
        
        try:
            allpair_file_path: str = [
                        file for file in allpair_list
                        if self.fix_chr(self.chr_num.strip(".")) in file and self.identifier in file
                    ][0]

        except IndexError:

            print(
                    f"There was no allpair.txt file found for the variant, {self.identifier}"
                )
            # need to catch an error here
            return "None"

        return allpair_file_path
    
    @staticmethod
    def load_allpair_file(allpair_file_path: str) -> pd.DataFrame:
        """This function will load the allpair files into a dataframe
        Parameters
        __________
        allpair_file_path : str
            filepath to the allpair.txt file
            
        Returns
        _______
        pd.DataFrame
            returns a dataframe that has the pair1 value, the pair2 
            value, and then the carrier status for pair2
        """

        # Reading the file into a pandas dataframe
        return pd.read_csv(
            allpair_file_path,
            sep="\t",
            usecols=["pair_1", "pair_2", "carrier_status"])

    def filter_for_carriers(self, allpair_df: pd.DataFrame) -> pd.DataFrame:
        """Function to filter tha allpair_df for only where the carrier status = 1
        Parameters
        __________
        allpair_df : pd.DataFrame
            Dataframe that has the pair1 and pair2 and the carrier status
        
        Returns
        _______
        pd.DataFrame
            returns a dataframe where the pair 1 is in the carrier list
        """

        # making sure that second pair of grids are in the carrier list
        filtered_df: pd.DataFrame = allpair_df[allpair_df.carrier_status == 1]

        # making sure the first pair is in the carrier list
        filtered_df = filtered_df[filtered_df.pair_1.isin(self.iid_list)]

        return filtered_df

    def has_no_pairs(self, ind_in_networks_dict: Dict[str, List], output: str, analysis_type: str) -> Dict[str, List]:
        """Function to add the individuals who have no pairs to the output dictionary
        Parameters
        __________
        ind_in_networks_dict : Dict[str, List]
            dictionary that will be used to keep track of the variant_id, what 
            percentage of individuals are in the network, what number of carriers 
            are genotyped, and how many carriers are confirmed carriers.
        
        output : str
            string that list the filepath to the write the output 
            files to
        
        Returns 
        Dict[str, List]
            returns a dictionary containing information about the 
            variant and how many iids are confirmed carriers and 
            what percentage of them are in networks
        """
        if self.analysis_type == "phenotype":
            ind_in_networks_dict["variant_id"].append("N/A")
            ind_in_networks_dict["gene_name"].append(self.identifier)
            ind_in_networks_dict["percent_in_network"].append(0.00)
            ind_in_networks_dict["genotyped_carriers_count"].append(0)
            ind_in_networks_dict["confirmed_carrier_count"].append(0)
        else:
            ind_in_networks_dict["variant_id"].append(self.identifier)
            ind_in_networks_dict["gene_name"].append("N/A")
            ind_in_networks_dict["percent_in_network"].append(0.00)
            ind_in_networks_dict["genotyped_carriers_count"].append(len(self.iid_list))
            ind_in_networks_dict["confirmed_carrier_count"].append(0)

        # creating a dictionary to put information into
        carriers_in_network_dict: Dict[str, List] = {
            "IID": [],
            "In Network": [],
            "Network ID": [],
            "gene_name":[],
            "variant_id":[],
            "chr_num":[]
        }

        # iterating through each iid in the iid_list so that
        # we can record the iids that have no pairs and therefore 
        # are not in a network
        if analysis_type == "phenotype":
            for iid in self.iid_list:
                carriers_in_network_dict["IID"].append(iid)
                carriers_in_network_dict["In Network"].append(0)
                carriers_in_network_dict["Network ID"].append("N/A")
                carriers_in_network_dict["gene_name"].append(self.identifier)
                carriers_in_network_dict["variant_id"].append("N/A")
                carriers_in_network_dict["chr_num"].append(self.fix_chr(self.chr_num)[-2:])
        else:
            for iid in self.iid_list:
                carriers_in_network_dict["IID"].append(iid)
                carriers_in_network_dict["In Network"].append(0)
                carriers_in_network_dict["Network ID"].append("N/A")
                carriers_in_network_dict["gene_name"].append("N/A")
                carriers_in_network_dict["variant_id"].append(self.identifier)
                carriers_in_network_dict["chr_num"].append(self.fix_chr(self.chr_num)[-2:])

        # converting the above dictionary to a dataframe 
        self.network_carriers = pd.DataFrame.from_dict(
            carriers_in_network_dict)
        print("self.network_carriers")
        print(self.network_carriers)
        # need to write these variants to a dataframe
        # if the file already exist then need to append the new 
        # information without adding a header
        if os.path.exists(os.path.join(
            output, "network_groups.txt")) and os.stat(os.path.join(
                output, "network_groups.txt")) != 0:

            self.network_carriers.to_csv(os.path.join(
                output, "network_groups.txt"), 
                                        sep="\t",
                                        index=False,
                                        mode="a",
                                        header=None)
        
        # if the file doesn't exist already then you have to append
        # the new information with a header
        else:
            self.network_carriers.to_csv(os.path.join(
                output, "network_groups.txt"),
                                        sep="\t",
                                        index=False,
                                        mode="a")

        return ind_in_networks_dict
    
    @staticmethod
    def drop_empty_rows(filtered_df: pd.DataFrame) -> pd.DataFrame:
        """Function to check if there are any nul values in the filtered allpair_df and if there are then those rows get dropped
        Parameters
        __________
        filtered_df : pd.DataFrame
            dataframe that has been been filtered for carriers
        
        Returns
        _______
        pd.DataFrame
            dataframe that has no null values"""
        # This just drops any empty rows in the dataframe
        if filtered_df["pair_1"].isnull().any(
        ) or filtered_df["pair_2"].isnull().any():
            nan_value = float("NaN")
            filtered_df.replace("", nan_value, inplace=True)
            filtered_df.dropna(subset=["Pairs"], inplace=True)

        return filtered_df
                

    def carriers_in_network(self, subset_df: pd.DataFrame,
                            ind_in_networks_dict: dict) -> dict:
        '''This function tells the percent of carriers who are in these networks'''

        id1_set: set = set(subset_df.pair_1.values.tolist())

        id2_set: set = set(subset_df.pair_2.values.tolist())

        total_carriers_set: set = id1_set | id2_set

        carriers_in_network: int = sum(item in total_carriers_set
                                       for item in set(self.iid_list))


        percent_in_networks: float = round(carriers_in_network / len(self.iid_list),
                                           3)

        if self.analysis_type == "phenotype":
            ind_in_networks_dict["variant_id"].append("N/A")
            ind_in_networks_dict["gene_name"].append(self.identifier)
            ind_in_networks_dict["percent_in_network"].append(percent_in_networks)
            ind_in_networks_dict["genotyped_carriers_count"].append(0)
            ind_in_networks_dict["confirmed_carrier_count"].append(
                carriers_in_network)
        else:
            # adding values into the above dictionary
            ind_in_networks_dict["variant_id"].append(self.identifier)
            ind_in_networks_dict["gene_name"].append("N/A")
            ind_in_networks_dict["percent_in_network"].append(percent_in_networks)
            ind_in_networks_dict["genotyped_carriers_count"].append(len(self.iid_list))
            ind_in_networks_dict["confirmed_carrier_count"].append(
                carriers_in_network)
        
        # Need to make it so that there is a dataframe with the IID and then a 1 or 0 if it is in a network or not
        # figure out which carriers from the iid_list are in networks
        carriers_in_network_dict = {
            "IID": [],
            "In Network": [],
            "Network ID": []
        }

        # This for loop checks to see if the iid in the iid_list shares a segemnt. It it will then return a 1 if
        # the iid is a carrier or a 0 if it is not avaliable
        
        for iid in self.iid_list:

            bool_int = int(iid in total_carriers_set)

            carriers_in_network_dict["IID"].append(iid)
            carriers_in_network_dict["In Network"].append(bool_int)
            carriers_in_network_dict["Network ID"].append("N/A")

        # convert dictionary to Dataframe

        self.network_carriers = pd.DataFrame(
            carriers_in_network_dict,
            columns=["IID", "In Network", "Network ID"])
        
        return ind_in_networks_dict

    def draw_networks(self, reformated_df: pd.DataFrame, pairs_df: pd.DataFrame, output_path: str) -> tuple:
        """This function actually draws the networks. It takes the reformated dataframe from the isolate_ids functions. It will return the path to the output file and a dataframe of pairs"""

        ##########################################################

        # getting a list of all unique IDs in the Pair_id1 column to iterate through
        id_list = reformated_df['pair_1'].unique().tolist()

        # keeping a set of nodes that have already been visited
        nodes_visited_set = set()
        # Starting the network number at 1
        network_number = 1
        # Creating the nodes
        # iterating through each id1
        for id1 in id_list:

            # Limit the passed dataframe to only those that have a relationship with id1. Only need to Pair_id1 and
            # Pair_id2 columns to make nodes
            nodes_constructed = set()

            if id1 not in nodes_visited_set:

                segments_file_subset = reformated_df[
                    (reformated_df['pair_1'] == id1) |
                    (reformated_df['pair_2'] == id1)][["pair_1", "pair_2"]]

                ids1_set = set(segments_file_subset.pair_1.values.tolist())

                ids2_set = set(segments_file_subset.pair_2.values.tolist())

                ids_list = ids1_set | ids2_set

                # insert recursive function here
                self.generate_pair_list(ids_list, reformated_df, id1)

                # Subsetting the original dataframe for this full list
                full_subset_df = reformated_df[
                    (reformated_df.pair_1.isin(self.iid_list)) |
                    (reformated_df.pair_2.isin(self.iid_list))][[
                        "pair_1", "pair_2"
                    ]]

                # iterating through each row
                for row in full_subset_df.itertuples():

                    # breaking up the tuple into each pair
                    Pair_id1 = row[1]
                    Pair_id2 = row[2]
                    # This if statements only makes Pair_id1 a node if it is not already constructed
                    if Pair_id1 not in nodes_constructed:

                        # Keeping track of what nodes have been made
                        nodes_visited_set.add(Pair_id1)

                        nodes_constructed.add(Pair_id1)

                    # This if statement only makes the second node of Pair_id2 if it has not been made yet
                    if Pair_id2 not in nodes_visited_set:


                        # Keeps track of what nodes have been made
                        nodes_visited_set.add(Pair_id2)

                        nodes_constructed.add(Pair_id2)


                ###########################################################################
                # Next section is responsible for updating which IID is in which network
                # Updating the self.network_carriers dataframe so that the Network id gets set to a number

                    idx = self.network_carriers.index[
                        self.network_carriers["IID"] == Pair_id1]

                    self.network_carriers.loc[idx, ["Network ID"]] = str(
                        network_number)

                    idx = self.network_carriers.index[
                        self.network_carriers["IID"] == Pair_id2]

                    self.network_carriers.loc[idx, ["Network ID"]] = str(
                        network_number)

                ###########################################################################

                # This updates which network is being created
                network_number += 1
                # Adding values to the pair dictionary


        # Adding a column for the variant number and the chr number so that we can output just one file
        self.network_carriers = self.add_columns(self.network_carriers)
        print(self.network_carriers)
        # allpairs_df = self.add_columns(allpairs_df, variant, chr_num)

        # Writing the dataframe to a csv file
        if os.path.exists(os.path.join(
            output_path, "network_groups.txt")) and os.stat(os.path.join(
                output_path, "network_groups.txt")) != 0:

            self.network_carriers.to_csv(os.path.join(
                output_path, "network_groups.txt"),
                                        sep="\t",
                                        index=False,
                                        mode="a+",
                                        header=None)
        else:

            self.network_carriers.to_csv(os.path.join(
                output_path, "network_groups.txt"),
                                        sep="\t",
                                        index=False,
                                        mode="a")
        # the function returns the full path for the network_carriers dataframe and it also returnst eh allpairs_df
        return os.path.join(output_path, "network_groups.txt"), pairs_df

    def add_columns(self, dataframe: pd.DataFrame,
                    ) -> pd.DataFrame:
        '''This function will take the carrier dataframe and add two columns for the variant id and the chr_num. It will then return the dataframe'''

        
        if self.analysis_type == "phenotype":
            dataframe["gene_name"] = self.identifier
            dataframe["variant_id"] = "N/A"
        else:
            # adding a column for the variant id
            dataframe["gene_name"] = "N/A"
            dataframe["variant_id"] = self.identifier

        # adding a column for the chromosome number

        dataframe["chr_num"] = self.fix_chr(self.chr_num)[-2:]

        return dataframe
    
    def generate_pair_list(self, current_id_set: set, segment_df: pd.DataFrame,
                           current_id1: str):
        '''This function identifies all the potential pair ids for '''

        # First subset the dataframe for all the current ids
        new_df_subset:pd.DataFrame = segment_df[
            (segment_df['pair_1'].isin(current_id_set)) |
            (segment_df['pair_2'].isin(current_id_set))][["pair_1", "pair_2"]]

        # These next three lines use sets to avoid duplicate values
        # This gets all ids from the pair one column
        id1_set: set = set(new_df_subset.pair_1.values.tolist())
        # This gets all ids from the Pair two column
        id2_set: set  = set(new_df_subset.pair_2.values.tolist())
        # this line takes teh union of each set to get all unique values
        total_id_set = id1_set | id2_set

        # If the total_id_set equals the current_id_set then the function stops
        if total_id_set == current_id_set:
           
            # returns the set of ids
            self.iid_list: set = total_id_set
            
            return
        # If the two sets are not equal then it recursively calls the function again
        else:

            self.generate_pair_list(total_id_set, segment_df, current_id1)
