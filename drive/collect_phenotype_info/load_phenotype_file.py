import pandas as pd
import sys

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
    # trying to load in the phenotype file
    try:
        pheno_df: pd.DataFrame = pd.read_csv(pheno_filepath, sep="\t")
    except FileNotFoundError:
        print(f"The provided input file was not found at: {pheno_filepath}")
        sys.exit(1)
    # trying to load in the carriers file
    try:
        carriers_df: pd.DataFrame = pd.read_excel(carriers_file)
    except FileNotFoundError:
        print(f"The provided input file was not found at: {carriers_file}")
        sys.exit(1)

    return pheno_df, carriers_df

        