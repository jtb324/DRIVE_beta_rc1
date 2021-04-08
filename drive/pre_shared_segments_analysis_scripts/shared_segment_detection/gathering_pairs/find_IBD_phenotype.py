# This file will be a way to form the small.txt.gz file for all pairs within the provided ibd file for a certain range
import pandas as pd
import sys
import os
import multiprocessing as mp
from functools import partial

from .collect_shared_segments import generate_parameters, build_unique_id_dict, create_ibd_arrays, gather_pairs
import utility_scripts

# making a custom exception for when the iid_list is empty

class EmptyList(Exception):
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors
        print("Printing Errors: ")
        
def get_carriers(phenotype_carriers_df: pd.DataFrame) -> list:
    """Function to get the list of iids identified as having the phenotype
    Parameters
    __________
    phenotype_carriers_df : pd.DataFrame
        dataframe that has one column, IID, and list the grid iids
        for grids that were identified as carrying the phenotype of 
        interest

    Returns
    _______
    list
        returns a list of these carriers
    """
    try:
        iid_list: list = phenotype_carriers_df.IID.values.tolist()
    except KeyError:
        raise KeyError("Column IID not found in the phenotype carriers file.")

    if len(iid_list) == 0:
        raise EmptyList("Expected at least grid to be in the provided phenotype carriers file","len(iid_list) == 0")

    return iid_list

def create_dictionary(row: pd.Series, gene_dict: dict, ):
    """Function to create the dictionary for each gene"""

    gene_dict[row[0]] = {
        "chr": row[1],
        "start": row[2],
        "end": row[3],
        }

def gather_gene_info(gmap_df: str) -> dict:
    """Function that will return a dictionary of all the important info for the gene
    Parameters
    __________
    gmap_file : str
        Tab separated file containing information about the gene such as gene name, chromosome, start, and end position
    
    Returns
    _______
    dict
        returns a dictionary where the genes are the key and the 
        values are the chromosome number, the start position, and 
        end position
    """

    gene_dict: dict = {}

    gmap_df.apply(lambda row: create_dictionary(row, gene_dict), axis=1)

    return gene_dict

def get_ibd_file(ibd_list: list, chr_num) -> str:
    """Function to get the ibd_file from the ibd_list
    Parameters
    __________
    ibd_list : list
        a list containing all of the ibd_files from the specified directory
        
    chr_num : str
        string containing the chromosome number. This will be just a digit as a 
        string
    
    Returns
    _______
    str
        string that list the pathway to the ibd_file for the specified chromosome 
    """

    chr_num = "".join(["chr", str(chr_num),"."])

    return [file for file in ibd_list if chr_num in file][0]

def collect_IBD_segments(carrier_list: list, ibd_program:str, min_CM: str, ibd_file_list: list, output_path: str, gene_dict: dict, que_object, key: str):

    gene_info: dict = gene_dict[key]
    
    ibd_file: str = get_ibd_file(ibd_file_list, gene_info["chr"])

    # getting the indices of all the values of interest for the     
    # ibd_program used
    parameter_dict: dict = generate_parameters(ibd_program)

    # building the uniqID dict
    uniqID: dict = build_unique_id_dict(carrier_list)

    IBDdata, IBDindex = create_ibd_arrays()

    chr_num: str = gather_pairs(IBDdata, IBDindex, parameter_dict, ibd_file, uniqID, min_CM, que_object, output_path, ibd_program, gene_start=gene_info["start"], gene_end=gene_info["end"], gene_name=key) 

def run_parallel(gene_info_dict: dict, ibd_file_list: list,THREADS: int, min_CM: str, ibd_program: str, output: str, carrier_list: list):
    """function to run through the genes in parallel"""

    manager = mp.Manager()

    que = manager.Queue()

    pool = mp.Pool(int(THREADS))
    header:str = "gene\n"

    watcher = pool.apply_async(
            utility_scripts.listener,
            (que, "".join([output, "gene_target_failed.txt"]), header))

    func = partial(collect_IBD_segments, carrier_list, ibd_program, min_CM, ibd_file_list, output, gene_info_dict, que)

    pool.map(func, list(gene_info_dict.keys()))

    que.put("kill")

    pool.close()

    pool.join()


def gather_shared_segments(ibd_file_list: list, pheno_gmap_df:pd.DataFrame, phenotype_carriers_df: pd.DataFrame, output_path: str, ibd_program: str, min_CM: str, ibd_suffix: str, THREADS):
    """Function to get the shared segments for each pair within a gene of interest
    Parameters
    __________
    ibd_file_path : string
        filepath to the directory of the ibd files from either 
        hapibd or ilash
    
    pheno_gmap_df: pd.DataFrame
        dataframe containing information about the gene of interest and the start and end point as well as the chromosome that the gene is on in the format of chrXX where X is a digit
        
    phenotype_carriers_df : pd.DataFrame
        dataframe containing a single column called IID and the 
        list each grid identified as having the phenotype of 
        interest
        
    output_path : str
        a string listing the output path to write the file to
    
    ibd_program : str
        This is the ibd program used to get the shared segment data. Should be either ilash or hapibd"""

    # checking to make sure the output directory subdirectory "collected_pairs" exists
    output_dir: str = utility_scripts.check_dir(output_path, "collected_pairs/")
    # getting a list of ibd_files
    ibd_file_list: list = utility_scripts.get_file_list(ibd_file_list, ibd_suffix)

    # getting a list of grids that have the phenotype of interest
    carrier_list: list = get_carriers(phenotype_carriers_df)

    # need to generate a dictionary of all chromosomes, with their start and end point
    gene_dict: dict = gather_gene_info(pheno_gmap_df)

    run_parallel(gene_dict, ibd_file_list, THREADS, min_CM, ibd_program, output_dir, carrier_list)

