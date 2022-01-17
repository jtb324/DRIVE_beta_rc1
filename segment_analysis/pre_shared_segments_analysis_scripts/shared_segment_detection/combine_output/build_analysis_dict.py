import pandas as pd
# This script has the functions that will build the 
# analysis_info_dict that has the keys for the analysis and either 
# the variant position or the gene start and end point
def get_variant_pos(map_file: str, identifier: str) -> tuple:
    """Function to get the variant position out of the map file
    Parameters
    __________
    map_file : str
        file path to the map file formed by initially running  PLINK
        
    identifier : str
        string that will list the variant id
    
    Returns
    _______
    tuple
        returns a tuplethat contains the variant position and then the variant id
    """
    # loading the map file into a dataframe
    map_df: pd.DataFrame = pd.read_csv(map_file, sep="\t", header=None)
    
    # getting the position of the variant out of the dataframe
    variant_position: str = map_df[map_df[1] == identifier[:-2]][3].values.tolist()[0]

    return variant_position, identifier


def get_gene_pos(identifier: str, pheno_gmap_df: pd.DataFrame) -> tuple:
    """Function to get the gene start and end position from the pheno_gmap_df
    Parameters
    __________
    identifier : str
        this is the string which is the gene name
    
    pheno_gmap_df : pd.DataFrame
        this is the dataframe that tells information about the gene 
        such as the gene name, chromosome the gene is on, and gene 
        start and end
    
    Returns
    _______
    tuple
        tuple of the gene start and end positions
    """
    start_pos: str = str(pheno_gmap_df[pheno_gmap_df[0] == identifier][2].values.tolist()[0])

    end_pos: str = str(pheno_gmap_df[pheno_gmap_df[0] == identifier][3].values.tolist()[0])

    return start_pos, end_pos

def get_analysis_files(analysis_type: str, identifier: str, map_file: str = None, pheno_gmap_df: pd.DataFrame=None) -> dict:
    """Function to return a dictionary with the analysis type and either the variant position or the gene start and end point
    Parameters
    __________
    analysis_type : str
        string that has the analysis type. This will either be phenotype, gene, or blank
        
    identifier : str
        string that is either the variant id or the gene name
    
    map_file : str
        This is the filepath to the map file if the map_file_list is a key in the gathered_file_dict. If not this value is none
    
    pheno_gmap_df : pd.DataFrame
        dataframe that has the gene of interest as well as information like the chromosome and the start and end position. This value will be none by default

    Returns
    _______
    dict
        returns a dictionary where the keys are the analysis type, variant position, or the start and end point for the gene
    """

    if analysis_type == "phenotype":

        gene_start, gene_end = get_gene_pos(identifier, pheno_gmap_df)

        return {
            "analysis_type": analysis_type,
            "gene_start": gene_start,
            "gene_end": gene_end
        }
    else:
        var_pos, _ = get_variant_pos(map_file, identifier)

        return {
            "analysis_type": analysis_type,
            "variant_pos": var_pos
            }