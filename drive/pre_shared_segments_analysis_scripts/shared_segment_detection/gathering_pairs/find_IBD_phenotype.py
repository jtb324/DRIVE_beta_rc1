# This file will be a way to form the small.txt.gz file for all pairs within the provided ibd file for a certain range
import pandas as pd
import sys
import os

from .collect_shared_segments import generate_parameters, build_unique_id_dict, create_ibd_arrays, gather_pairs
import utility_scripts
# making a custom exception for when the iid_list is empty

class EmptyList(Exception):
    def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors
        print("Printing Errors: ")
        print(errors
        )
        
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
    print(row)
    gene_dict[row[0]] = {
        "chr": row[1],
        "start": row[2],
        "end": row[3],
        }

def gather_gene_info(gmap_file: str) -> dict:
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
    # load the gmap_file into a dataframe
    gmap_df: pd.DataFrame = pd.read_csv(gmap_file, sep="\t")

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

    chr_num = "".join(["chr", chr_num,"."])

    return [file for file in ibd_list if chr_num in file][0]



def gather_shared_segments(ibd_file_list: list, pheno_gmap_df:pd.DataFrame, phenotype_carriers_df: pd.DataFrame, output_path: str, ibd_program: str, min_CM: str, ibd_suffix: str):
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
    
    # getting a list of ibd_files
    ibd_file_list: list = utility_scripts.get_file_list(ibd_file_list, ibd_suffix)

    # getting a list of grids that have the phenotype of interest
    carrier_list: list = get_carriers(phenotype_carriers_df)

    # need to generate a dictionary of all chromosomes, with their start and end point
    gene_dict: dict = gather_gene_info(pheno_gmap_df)

    for key in gene_dict:
        # pulling out the inner dictionaries from the gene key
        gene_info: dict = gene_dict[key]
    
        ibd_file: str = get_ibd_file(ibd_file_list, gene_info["chr"])
        # getting the indices of all the values of interest for the     
        # ibd_program used
        parameter_dict: dict = generate_parameters(ibd_program)

        # building the uniqID dict
        uniqID: dict = build_unique_id_dict(carrier_list)

        
        IBDdata, _ = create_ibd_arrays()

        chr_num: str = gather_pairs(IBDdata, IBDdata, parameter_dict, ibd_file, uniqID, min_CM) 

