import pandas as pd
import sys
import os


def get_var_pos(map_file: str, variant_id: str) -> int:
    '''This function find the position of the variant id and return that as a string'''

    # This reads in the columns from the map file that contains the variant probe id and the pos
    map_df: pd.DataFrame = pd.read_csv(map_file,
                                       sep="\t",
                                       header=None,
                                       usecols=[1, 3])

    # getting the value
    # print("int the get var pos function")
    # print(variant_id)
    # print(map_df)
    # print(type(map_df[map_df[1].isin([variant_id[:-2]])]))
    pos: int = int(map_df[map_df[1].isin([variant_id[:-2]])][3].values[0])

    return pos


def get_phase_number(pair_df_row: pd.DataFrame, program: str) -> list:
    """function to get the phase from each individual 
    Parameters
    __________
    pair_df_row : pd.DataFrame
        dataframe for the row that contains the pair of interest within the 
        ibd files

    program : str
        string that indicates which ibd program was used. This value will be 
        either ilash or hapibd

    Returns
    _______
    list
        returns a list where the first value is the phasing of the first 
        haplotype and the second values is the phasing of the second one
    """
    if program == "ilash":
        pairs: pd.DataFrame = pair_df_row[[1, 3]]

        pair1: str = pairs[1].values[0]
        pair2: str = pairs[3].values[0]

        phase1: str = pair1.split("_")[1]
        phase2: str = pair2.split("_")[1]

        return [phase1, phase2]

    elif program == "hapibd":
        pairs: pd.DataFrame = pair_df_row[[1, 3]]
        # print("in the get phaser number function")
        # print(pairs)
        phase1: str = pairs[1].values[0]
        phase2: str = pairs[3].values[0]

        return [phase1, phase2]
    else:
        print(
            "The provided program was not recognized. Please use output from either hapibd or ilash"
        )
        sys.exit(1)


def search_ibd_file(ibd_file: str, pair1: str, pair2: str, var_pos: int,
                    program: str):
    """function to find the start,end postions and segment length in ibd files for a pair for a variant
    Parameters
    __________
    ibd_file: str 
        string contain the file path to the hapibd or ilash output files

    pair1 : str
        string that contains the GRID id for the first pair

    pair2 : str
        string that contains the GRID id for the second pair

    var_position : int
        base position for the variant of interest based on the map file from plink
    
    program : str
        value will be either hapibd or ilash depend on which ibd_file is being used
    
    Return
    ______
    dict
        dictionary containing the start and end positions as well as the segment length
    """

    if program == "ilash":
        for chunk in pd.read_csv(ibd_file,
                                 sep="\t",
                                 chunksize=100000,
                                 header=None):

            filtered_chunk: pd.DataFrame = chunk[(chunk[0] == pair1)
                                                 & (chunk[2] == pair2)]

            if not filtered_chunk.empty:
                print("in the search ibd file function")
                print(filtered_chunk)
                print(pair1)
                print(pair2)
                start_pos: int = int(chunk[5].values[0])
                end_pos: int = int(chunk[6].values[0])

                if start_pos <= var_pos and end_pos >= var_pos:
                    length: str = str(chunk[9].values[0])

                    # get the two phase string
                    phase_list: list = get_phase_number(
                        filtered_chunk, "ilash")

                    return {
                        "start": start_pos,
                        "end": end_pos,
                        "length": length,
                        "phase1": phase_list[0],
                        "phase2": phase_list[1]
                    }
        if filtered_chunk.empty:
            return {
                "start": "N/A",
                "end": "N/A",
                "length": "N/A",
                "phase1": "N/A",
                "phase2": "N/A"
            }

    elif program == "hapibd":

        for chunk in pd.read_csv(ibd_file,
                                 sep="\t",
                                 chunksize=10000,
                                 header=None):

            filtered_chunk: pd.DataFrame = chunk[(chunk[0] == pair1)
                                                 & (chunk[2] == pair2)]

            if not filtered_chunk.empty:
                start_pos: int = int(chunk[5].values[0])
                end_pos: int = int(chunk[6].values[0])

                if start_pos <= var_pos and end_pos >= var_pos:
                    length: str = str(chunk[7].values[0])

                    # get the two phase string
                    phase_list: list = get_phase_number(
                        filtered_chunk, "hapibd")

                    return {
                        "start": start_pos,
                        "end": end_pos,
                        "length": length,
                        "phase1": phase_list[0],
                        "phase2": phase_list[1]
                    }
        # returns a dictionary of N/As if the pair was not found within the file
        if filtered_chunk.empty:
            return {
                "start": "N/A",
                "end": "N/A",
                "length": "N/A",
                "phase1": "N/A",
                "phase2": "N/A"
            }

    else:
        print(
            "The provided ibd program was not recognized. Please provide either a hapibd or ilash file"
        )
        sys.exit(1)


def get_segment_lens(pair1: str, pair2: str, map_file: str, variant_id: str,
                     hapibd_file: str, ilash_file: str):

    var_pos: int = get_var_pos(map_file, variant_id)

    hapibd_segment_dict: dict = search_ibd_file(hapibd_file, pair1, pair2,
                                                var_pos, "hapibd")
    ilash_segment_dict: dict = search_ibd_file(ilash_file, pair1, pair2,
                                               var_pos, "ilash")

    return hapibd_segment_dict, ilash_segment_dict
