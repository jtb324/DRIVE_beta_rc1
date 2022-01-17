import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import utility_scripts
from typing import List, Dict


def individual_distributions(networks_df: pd.DataFrame, analysis_type: str, output_dir: str) -> None:
    """Function that will create the graph for each gene/variant
    network_df : pd.DataFrame
        dataframe that results from loading the network_groups.txt file into memory

    analysis_type : str
        string that tells where the rare variant "gene" based analysis 
        is being used or the phenotype driven approach

    output_dir : str
        directory to output files into
    """
    # getting a list of genes/varaints
    identifier_handler: Dict = {
            "phenotype": list(set(networks_df["gene_name"].values.tolist())), "gene": list(set(networks_df["variant_id"].values.tolist()))
        }
    
    identifier_list: List[str] = identifier_handler[analysis_type]

    size_list: List[int] = []

    for id in identifier_list:
        subset_handler: Dict = {
                "phenotype": networks_df[networks_df["gene_name"].isin([id])], "gene": networks_df[networks_df["variant_id"].isin([id])]
            }
        
        network_subset: pd.DataFrame = subset_handler[analysis_type]

        # get the list of network sizes
        network_sizes: List[int] = list(network_subset["Network ID"].value_counts())

        # keeping track of all the distributions so a total
        # distribution can be determined
        for element in network_sizes:

            size_list.append(element)

        fig = plt.figure()

        plt.hist(network_sizes, edgecolor='black') 
        plt.title("".join(["Distribution of network sizes for ", id]))
        plt.ylabel("Count")
        plt.xlabel("Size of Network")
        plt.savefig(os.path.join(output_dir, "".join([id, "_network_size_distribution.png"])))

    fig = plt.figure()

    plt.hist(size_list, edgecolor='black') 
    plt.title("".join(["Distribution of network sizes for all genes/variants"]))
    plt.ylabel("Count")
    plt.xlabel("Size of Network")
    plt.savefig(os.path.join(output_dir, "".join(["Total_distributions.png"])))

# def total_distribution(networks_df: pd.DataFrame, output_dir: str) -> None:
#     """Function that will create the graph for each gene/variant
#     network_df : pd.DataFrame
#         dataframe that results from loading the network_groups.txt file into memory

#     output_dir : str
#         directory to output files into
#     """
#     print(networks_df["Network ID"].value_counts())
#     size_distribution_list: List[int] = list(networks_df["Network ID"].value_counts())
#     print(size_distribution_list)
#     print(max(size_distribution_list))
#     fig = plt.figure()

#     plt.hist(size_distribution_list, edgecolor='black') 
#     plt.title("".join(["Distribution of network sizes for all genes/variants"]))
#     plt.ylabel("Count")
#     plt.xlabel("Size of Network")
#     plt.savefig(os.path.join(output_dir, "".join(["Total_distributions.png"])))

def create_distributions(network_dir: str, analysis_type: str) -> None:
    """Function that will create the distribution of sizes for each network
    Parameters
    _________
    network_dir : str
        directory that has the network_groups.txt file

    analysis_type : str
        string that tells where the rare variant "gene" based analysis 
        is being used or the phenotype driven approach
    """
    print("Drawing plots of the distributions of the size of each network")
    # first create a new subdirectory to put the graphs into
    network_graphs_dir: str = utility_scripts.check_dir(network_dir, "network_distributions")

    networks_file: str = os.path.join(network_dir, "network_groups.txt")

    networks_df: pd.DataFrame = pd.read_csv(networks_file, sep="\t")

    # getting the individual distributions
    individual_distributions(networks_df, analysis_type, network_graphs_dir)

    # # getting the total distributions
    # total_distribution(networks_df, network_graphs_dir)

    
