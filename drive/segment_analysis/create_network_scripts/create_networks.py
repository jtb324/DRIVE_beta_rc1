
import logging
import re
import os
import pandas as pd
from typing import List, Dict

from segment_analysis.create_network_scripts.network_creator_class import Network_Maker, Network_Prep
import utility_scripts

def create_readme(output_path: str):
    '''This function will create a readme file for the specified directory'''
    readme = utility_scripts.Readme("_README.md", output_path)
    readme.rm_previous_file()
    readme.write_header("networks/")
    readme.create_date_info()
    readme.add_line(utility_scripts.networks_body_text)


def add_nopair_variants(ind_in_networks_dict: dict, iid_list: list,
                        variant: str) -> dict:
    '''This function will add the variants that have no pairs of confirmed 
    carriers to the list'''

    ind_in_networks_dict["variant_id"].append(variant)
    ind_in_networks_dict["percent_in_network"].append(0.00)
    ind_in_networks_dict["genotyped_carriers_count"].append(len(iid_list))
    ind_in_networks_dict["confirmed_carrier_count"].append(0)

    return ind_in_networks_dict


def create_networks(allpair_file_dir: str, networks_dir: str, analysis_type: str, confirmed_carriers_file: str):
    """main function that will be used to create networks of files"""
    create_readme(networks_dir)

    # checking to see if the next three files exist from a previous run and if they do then the progam removes them
    utility_scripts.check_file(os.path.join(networks_dir, "network_groups.csv"))

    utility_scripts.check_file(os.path.join(networks_dir, "/pairs_in_networks.csv"))

    utility_scripts.check_file(os.path.join(networks_dir, "missing_allpairs.txt"))

    # Creating an empty dataframe that the pairs will be written to at the end in the
    # draw networks function
    pairs_df: pd.DataFrame = pd.DataFrame()

    # Getting all of the allpair files into a list
    allpair_file_list: list = utility_scripts.get_file_list(
        allpair_file_dir, "*.allpair.txt")

    # Need to refactor so that this doesn't rel on the variant_file_dir
    network_drawer: Network_Prep = Network_Prep(
        allpair_file_dir, networks_dir, analysis_type)
    
    network_drawer.determine_carriers(confirmed_carriers_file)
    
    
    # Getting a list of all the carrier files

    # Creating a dictionary that will be used to record useful information
    ind_in_networks_dict = {
        "variant_id": [],
        "gene_name": [],
        "percent_in_network": [],
        "genotyped_carriers_count": [],
        "confirmed_carrier_count": []
    }

    # creating a list that will be used to keep track of which variants have no pairs
    no_pair_carry_var_list: list = []


        # TODO: At this point the network_drawer has an attribute iid_dict that has the carriers for each variant/gene in for the chromosome 
        # iterating through each variant and getting the list of carriers
    for chromo_num in network_drawer.iid_dict.keys():
        
        inner_dict: Dict[str, List[str]] = network_drawer.iid_dict[chromo_num]
        # creating an object that has the chromosome number, the variant id/gene 
        # name, and the iid list as attributes

        for identifier in inner_dict: 
            # pulling out the list of iids that are associated with the different 
            # genes/variant ids
            iid_list: List[str] = inner_dict[identifier]

            # creating an object that will determine who is in a network
            network_maker: Network_Maker = Network_Maker(chromo_num, identifier, iid_list, analysis_type)

            allpair_file: str = network_maker.find_allpair_file(allpair_file_list)

            if allpair_file == "None":

                # TODO: add a function to write these variants to a file so that you can identify them

                continue
            
            # Loading the dataframe
            allpair_df: pd.DataFrame = network_maker.load_allpair_file(
                allpair_file)

            # filtering the allpair_df for only rows were both individuals are carriers
            filtered_allpair_df: pd.DataFrame = network_maker.filter_for_carriers(
                allpair_df)

            # if the dataframe is empty then the loop needs to skip that variant
            if filtered_allpair_df.empty:
                print(
                    f"There were no pairs found where both individuals were carriers for the variant {network_maker.identifier}"
                )

                # appending the variant/gene that has no pairs to the 
                # no_pair_carry_var_list
                no_pair_carry_var_list.append(network_maker.identifier)

                # adding the variants/gene that has no pairs to the ind_in_networks_dict
                ind_in_networks_dict = network_maker.has_no_pairs(ind_in_networks_dict, networks_dir)
                # prints the message and then moves onto the next iteration
                continue

            filtered_allpair_df = network_maker.drop_empty_rows(filtered_allpair_df)
            # This class method will determine the percentage of carriers in each network for each variant
            ind_in_networks_dict = network_maker.carriers_in_network(
                filtered_allpair_df, ind_in_networks_dict
                )

            groups_output_path, pairs_df = network_maker.draw_networks(
                filtered_allpair_df, pairs_df, networks_dir)

    ind_in_network_df: pd.DataFrame = pd.DataFrame.from_dict(
        ind_in_networks_dict)

    ind_in_network_df.to_csv(os.path.join(
        networks_dir, "percent_confirmed_carriers.csv"),
                             sep=",",
                             index=False)
