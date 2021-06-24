import pandas as pd

def filter_to_individual_in_uniqID(df_chunk: pd.DataFrame, uniqID: dict, id1_indx: int, id2_indx: int) -> pd.DataFrame:
    """Function to filter the dataframe chunk to only those rows where the id1 or id2 is 
    in the uniqID dictionary
    Parameters
    __________
    df_chunk : pd.DataFrame
        dataframe that is a chunk of the total ibd shared segment file from either 
        ilash or hapibd
    
    uniqID : dict 
        dictionary contain keys of iids that are indicated as carriers by the MEGA
        Ex array
    
    id1_indx : int
        integer tell the index that id1 is within the df_chunk

    id2_indx : int
        integer tell the index that id2 is within the df_chunk
    Returns
    _______
    pd.DataFrame
        returns a dataframe where either id1 or iƒd2 are in the uniqID dictionary
    """

    return df_chunk[(df_chunk[id1_indx].isin(uniqID)) |
                                    (df_chunk[id2_indx].isin(uniqID))]

def filter_to_greater_than_3_cm(df_chunk: pd.DataFrame, cM_indx: int, min_CM: int, 
                        unit_indx: int = None) -> pd.DataFrame:
    """Function to filter the df_chunk for values where the shared segment length is 
    greater than 3 cM 
    Parameters
    __________
    df_chunk : pd.DataFrame
        dataframe chunk from the ibd file where the id1 and id2 are in uniqID dict

    cM_indx : int
        integer that tells which column the cM value is in the df_chunk

    min_CM : int
        integer that tells the minimum centimorgan threshold 
    
    unit_indx : int
        this value tells what column the unit can be found. This value is only valid 
        for GERMLINE and will by default be none
    
    Returns
    _______
    pd.DataFrame
        returns a dataframe that is filtered for shared_segments of >= the minimum 
        centimorgan threshold
    """
    # This is reducing the dataframe to only pairs greater than min_cM threshold

    chunk_greater_than_3_cm = df_chunk[(
        df_chunk[cM_indx] >= min_CM)]
    # Need to check if the unit doesn't equal cM. This only applies in the case of germline
    if unit_indx:
        chunk_greater_than_3_cm = chunk_greater_than_3_cm[
            chunk_greater_than_3_cm[unit_indx] == "cM"]

    return chunk_greater_than_3_cm

def filter_for_correct_base_pair(df_chunk: pd.DataFrame, str_indx: int, end_indx: int, base_pair: int):
    """Function to filter chunk for values where the base pair is within the start and 
    end point
    Parameters
    __________
    df_chunk : pd.DataFrame
        dataframe that was filter by the filter_to_greater_than_3_cm
    
    str_indx : int
        integer that tells what column the start position is in the df_chunk
    
    end_indx : int
        integer that tells what column the end position is in the df_chunk
    
    base_pair : int
        integer that tells what the base pair position for the variant is

    Returns
    _______
    pd.DataFrame
        returns the filtered dataframe 
    """

    return df_chunk[(df_chunk[str_indx] < base_pair)
            & (df_chunk[end_indx] > base_pair)]

def filter_for_gene_site(df_chunk: pd.DataFrame, ibd_str_indx: int, ibd_end_indx: int, gene_start: int, gene_end: int) -> pd.DataFrame:
    """Function to filter tha dataframe chunk for sites where the 
    shared segment is partially in the gene or the the entire gene 
    is in the shared segment
    Parameters
    __________
    df_chunk : pd.DataFrame
        dataframe containing a chunk of information from the 
        ibd_file. This will be the output from the 
        filter_to_greater_than_3_cm function
    
    ibd_str_indx : int
        integer that gives the index of the ibd start position

    ibd_end_indx : int
        integer that give index of the ibd end position shared 
        segment
    
    gene_start : int
        integer that tells the base position of where the gene begins
    
    gene_end : int
        integer that tells the base position of where the gene ends
    
    Returns
    _______
    pd.DataFrame
        returns a pandas dataframe that is filtered for sites where part of the shared segment is within the gene or the gene is within the shared segment
    """

    
    # checking for the case where the shared segment end is before the gene end or where the shared segment start is after the gene start or where the gene is within the shared segment

    filtered_chunk: pd.DataFrame = df_chunk[((df_chunk[ibd_end_indx].astype(int) <= gene_end) & (df_chunk[ibd_end_indx].astype(int) >= gene_start))| ((df_chunk[ibd_str_indx].astype(int) >= gene_start) & (df_chunk[ibd_str_indx].astype(int) <= gene_end)) | ((df_chunk[ibd_str_indx].astype(int) <= gene_start) & (df_chunk[ibd_end_indx].astype(int) >= gene_end))]
    
    return filtered_chunk

def filter_for_matches(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Function that will filter the dataframe to make sure the pair1 id and pair 2 id are not the same
    Parameters
    __________
    dataframe : pd.DataFrame
        dataframe to be filtered
    
    Returns
    _______
    pd.DataFrame
        returns a dataframe that has been filtered so that each iids are different
    """

    return dataframe[dataframe[0] != dataframe[2]]

