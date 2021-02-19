# importing modules
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
import sys
from functools import partial
import shutil

# importing module from another file
import haplotype_segments_analysis
import utility_scripts


def get_variant_list(haplotype_len_df: pd.DataFrame) -> list:
    """
    Parameters
    __________
    haplotype_len_df : pd.DataFrame
        This parameter is a dataframe that contains information about the
        start and endpoints of each haplotype as well as the total length
        of the segment for both ilash and  hapibd

    Returns
    _______
    List
        Returns a list containing the different variant id within the
        dataframe
    """

    return list(set(haplotype_len_df.variant_id.values.tolist()))


def get_network_ids(haplotype_subset_df: pd.DataFrame,
                    variants_id: str) -> list:
    """
    Parameters
    __________
    haplotype_subset_df : pd.DataFrame
        pandas dataframe containing information about the start and endpoints
        of each haplotype for a specific variant as well as the toal length of 
        the segment for both ilash and hapibd and

    variant_id : str
        string containing the variant id

    Returns
    _______
    list
        returns a list of network ids. This list should only contains integers

    """

    return list(set(haplotype_subset_df.network_id.values.tolist()))


def get_chr_num(haplotype_info_df: pd.DataFrame) -> str:
    """
    Parameters
    __________
    haplotype_info_df : pd.DataFrame
        dataframe containing information about the start and endpoints
        for each pair's haplotype as well as the total length of the
        shared segment for both hapibd and ilash

    Returns
    _______
    str
        returns a string giving the chromosome number
    """
    chr_str: str = haplotype_info_df.chr.values.tolist()[0]

    return chr_str


def get_start_and_end(df_subset: pd.DataFrame, start_column_name: str,
                      end_column_name: str) -> tuple:
    """
    Parameters
    __________
    df_subset : pd.DataFrame
        Dataframe containing haplotype start and endpoint information for pairs
        in a specific network for a specific variant for both

    start_column_name : str
        String that gives the column name for the the column containing start
        positions of each haplotype. This value will be either hapibd_start or
        ilash start

    end_column_name : str
        String that gives the column name for the the column containing end
        positions of each haplotype. This value will be either hapibd_end or
        ilash end

    Returns
    _______
    tuple
        tuple containing two lists where one contains the start positions and
        another contains the end positions
    """
    start_positions_list: list = df_subset[start_column_name]

    start_positions_list = start_positions_list.dropna().values.tolist()

    end_positions_list: list = df_subset[end_column_name]

    end_positions_list = end_positions_list.dropna().values.tolist()
    return start_positions_list, end_positions_list


def get_median_value(positions_list: list) -> int:
    """
    Parameters
    __________
    positions_list : list
        List contain either the start or end positions of the shared segments
        for each pair in an network

    Returns
    int
        returns an integer of the median value of the list
    """
    positions_list.sort()

    midpoint: int = len(positions_list) // 2

    median = (positions_list[midpoint] + positions_list[~midpoint]) / 2

    return int(median)


def determine_haplotype_start_and_end(haplotype_len_df: pd.DataFrame,
                                      positions_dict: dict, variant_id: str,
                                      network_id: str):
    """
    Parameters
    __________
    haplotype_len_df : pd.DataFrame
        pandas dataframe containing information about the start and endpoints
        of each haplotype as well as the toal length of the segment for both
        ilash and hapibd and

    positions_dict : dict
        Dictionary that will be filled with the start and end positions for 
        each variant

    variant_id : str
        The id of the variant of interest

    network_id : str
        The id for the specific network
    """
    # subsetting the dataframe down to just the specific network
    haplotype_network_subset_df: pd.DataFrame = haplotype_len_df[
        haplotype_len_df.network_id == network_id]

    hapibd_start_list, hapibd_end_list = get_start_and_end(
        haplotype_network_subset_df, "hapibd_start", "hapibd_end")

    ilash_start_list, ilash_end_list = get_start_and_end(
        haplotype_network_subset_df, "ilash_start", "ilash_end")

    # getting the median value for hapibd and ilash start and end lists
    if hapibd_start_list and hapibd_end_list:
        hapibd_start_median: int = get_median_value(hapibd_start_list)
        hapibd_end_median: int = get_median_value(hapibd_end_list)

        positions_dict[variant_id][network_id]["hapibd"] = {
            "start": hapibd_start_median,
            "end": hapibd_end_median
        }

    if ilash_start_list and ilash_end_list:
        ilash_start_median: int = get_median_value(ilash_start_list)
        ilash_end_median: int = get_median_value(ilash_end_list)

        positions_dict[variant_id][network_id]["ilash"] = {
            "start": ilash_start_median,
            "end": ilash_end_median
        }
    # adding values for the hapibd and ilash start and endpoints into the
    # positions_dict


def form_dict(variant_list: list, haplotype_len_df: pd.DataFrame) -> dict:
    """
    Parameters
    __________
    variant_list : list
        list containing the id of every variant in the haplotype info files

    haplotype_len_df : pd.DataFrame
        pandas dataframe containing information about the start and endpoints
        of each haplotype as well as the toal length of the segment for both
        ilash and hapibd and

    Returns
    _______
    dict
        returns a nested dictionary that contains the median start and 
        endpoints for each network in each variant 
    """

    positions_dict: dict = dict()

    # TODO: Need to create a listener function that watches the cue and writes to
    # file when something is added to it
    # iterate through each variant
    for variant in variant_list:

        # subsetting the haplotype_len_df for a specific variant
        haplotype_subset_df: pd.DataFrame = haplotype_len_df[
            haplotype_len_df.variant_id == variant]

        # get the different networks within the dataframe
        networks_list: list = get_network_ids(haplotype_subset_df, variant)

        chr_num: str = get_chr_num(haplotype_subset_df)

        positions_dict[variant] = {
            id: {
                "chr": chr_num,
                "hapibd": None,
                "ilash": None
            }
            for id in networks_list
        }

        for network_id in networks_list:
            determine_haplotype_start_and_end(haplotype_subset_df,
                                              positions_dict, variant,
                                              network_id)

    return positions_dict


def rm_files(file_path: str):

    # get rid of the .ped ending first
    general_file_path: str = file_path[:-4]

    # removing the files iteratively
    for ending in [".ped", ".log", ".map"]:

        os.remove("".join([general_file_path, ending]))


def write_to_file(haplotype_str: str, output_path: str, variant: str,
                  chr_num: str, network_id: str, ibd_program: str,
                  position_list: list):
    """Function to write the haplotype string to a file 
    Parameters
    __________
    haplotype_str : str 
        string contain the most probable allele at each position

    output_path : str 
        filepath to write the output to 

    variant : str 
        The variant id that the haplotype is associated with 

    chr_num : str
        chromsome number that the haplotype is located on 

    network_id : str 
        The id of the specific network 

    position_list : list
        list of the start and end position of the haplotype

    """
    with open(output_path, "a+") as output_file:
        if os.path.getsize(output_path) == 0:
            output_file.write(
                "variant_id\tchr\tnetwork\tprogram\tstart_pos\tend_pos\thaplotype\n"
            )

        output_file.write(
            f"{variant}\t{chr_num}\t{network_id}\t{ibd_program}\t{position_list[0]}\t{position_list[1]}\t{haplotype_str}\n"
        )


def gather_haplotypes(haplotype_len_filepath: str, output: str,
                      binary_file: str, pop_info: str, pop_code: str) -> str:
    """
    Parameters
    ----------
    haplotype_len_filepath : str
        This parameter list the filepath to the file that contains the lengths 
        of each haplotype for each pair

    output : str 
        string listing the filepath to the output directory for the function

    binary_file : str
        string listing the directory to the bed and bim files used by PLINK

    pop_info : str
        string to the directory containing a file that list the population
        demographic for grids IIDs. This demographic codes should be based on 
        1000 genomes definition

    pop_code : str
        string listing the demographic code that the user wants to filter for

    Returns 
    -------
    str
        Returns a string containing the probable haplotype

    """

    # Loading the file containing haplotype lengths into a pandas dataframe
    haplotype_len_df: pd.DataFrame = pd.read_csv(haplotype_len_filepath, "\t")

    # generating a list of all the variants within the dataframe
    variants_list: list = get_variant_list(haplotype_len_df)

    start_end_dict: dict = form_dict(variants_list, haplotype_len_df)

    # creating a directory to put the plink files into
    try:
        os.mkdir("".join([output, "temp_plink_files/"]))
    except FileExistsError:
        pass

    plink_output_path: str = "".join([output, "temp_plink_files/"])

    output_path: str = "".join([output, "haplotypes/"])

    try:
        os.mkdir(output_path)
    except FileExistsError:
        pass

    output_file_path: str = "".join([output_path, "network_haplotypes.txt"])

    for variant, inner_dict in start_end_dict.items():

        for network_id, segment_dicts in inner_dict.items():

            chr_num: str = str(segment_dicts["chr"])
            # getting the hapibd/ilash start and endpoints from the dictionary
            if segment_dicts.get("hapibd"):
                hapibd_start: int = segment_dicts["hapibd"].get("start")
                hapibd_end: int = segment_dicts["hapibd"].get("end")

                hapibd_ped_file_path: str = haplotype_segments_analysis.get_plink_haplotype_str(
                    binary_file, str(hapibd_start), str(hapibd_end),
                    plink_output_path, "hapibd", chr_num, variant,
                    str(network_id))

                hapibd_haplotype_str: str = haplotype_segments_analysis.compare_haplotypes(
                    hapibd_ped_file_path, pop_info, pop_code)

                write_to_file(hapibd_haplotype_str, output_file_path, variant,
                              chr_num, network_id, "hapibd",
                              [hapibd_start, hapibd_end])

                rm_files(hapibd_ped_file_path)
                hapibd_ped_file_path = None

            if segment_dicts.get("ilash"):
                ilash_start: int = segment_dicts["ilash"].get("start")
                ilash_end: int = segment_dicts["ilash"].get("end")

                ilash_ped_file_path: str = haplotype_segments_analysis.get_plink_haplotype_str(
                    binary_file, str(ilash_start), str(ilash_end),
                    plink_output_path, "ilash", chr_num, variant,
                    str(network_id))

                ilash_haplotype_str: str = haplotype_segments_analysis.compare_haplotypes(
                    ilash_ped_file_path, pop_info, pop_code)

                write_to_file(ilash_haplotype_str, output_file_path, variant,
                              chr_num, network_id, "ilash",
                              [ilash_start, ilash_end])

                rm_files(ilash_ped_file_path)
                ilash_ped_file_path = None
                # TODO: need to add a file that catches if the haplotype string is not created for ilash or hapibd

                shutil.rmtree()
