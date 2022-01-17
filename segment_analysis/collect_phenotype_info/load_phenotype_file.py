import pandas as pd
import sys
import logging

def load_pheno_file(pheno_filepath: str, carriers_file: str) -> tuple:
    """Function to load the phenotype file into the program
    Parameters
    __________
    pheno_filepath : str 
        file path to the tab separated file that will load the phenotype file. This file will contain information about the chromosomes of interest and the position
    carriers_file : str 
        sfilepath to the input file that contains a list of grids that were identified as carrying a disease based off of Phenotype
    
    Returns
    _______
    tuple
        returns a tuple of two dataframes when the above files are loaded into a dataframe
    """
    # getting the main logger
    logger = logging.getLogger(__name__)

    # trying to load in the phenotype file
    try:

        pheno_df: pd.DataFrame = pd.read_csv(pheno_filepath, sep="\t", header=None)

    except FileNotFoundError:
        print(f"The provided input file was not found at: {pheno_filepath}")

        logger.error(f"The provided input file was not found at: {pheno_filepath}")
        sys.exit(1)
    # trying to load in the carriers file
    try:
        carriers_df: pd.DataFrame = pd.read_excel(carriers_file)
        
    except FileNotFoundError:
        print(f"The provided input file was not found at: {carriers_file}")
        logger.error(f"The provided input file was not found at: {carriers_file}")
        sys.exit(1)

    return pheno_df, carriers_df

def check_carrier_format(carriers_df: pd.DataFrame, logger) -> None:
    """Function to check if the provided file with carriers of the phenotype has two columns called IID and gene
    Parameters
    __________
    carriers_df : pd.DataFrame
        dataframe that list the carriers for a specific 
        phenotype and what gene to associate each 
        carrier with
    
    logger = logging.Logger
        logger object to record error messages
    """
    if "gene" not in carriers_df.columns or "IID" not in carriers_df.columns:
        print("The provided file of carriers does not have the correct columns. The file needs to have the columns 'IID' and the column 'gene'")
        logger.error("The provided file of carriers does not have the correct columns. The file needs to have the columns 'IID' and the column 'gene'")
    
        sys.exit(1)