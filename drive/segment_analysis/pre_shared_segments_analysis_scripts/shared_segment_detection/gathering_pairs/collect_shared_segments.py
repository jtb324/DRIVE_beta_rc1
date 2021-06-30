#!/usr/bin/python
import gzip
import multiprocessing as mp
import os
import utility_scripts
import pandas as pd

from typing import Union, Tuple, Dict

from ..generate_indx_dict.generate_dict import Germline_Indices, Ilash_Indices, Hapibd_Indices
from .filtering_functions import filter_to_greater_than_3_cm, filter_to_individual_in_uniqID, filter_for_correct_base_pair, filter_for_gene_site


####################################################################################################


class newPOS:
    __slots__ = 'add', 'rem'

    def __init__(self, add, rem):
        self.add = add
        self.rem = rem

def generate_parameters(ibd_program: str) -> dict:
    """Function to generate a dictionary of the indices for the parameters 
    that you have
    Parameters
    __________
    ibd_program : str
        string containing the ibd program that the input is coming from. 
        This value should be ilash, hapibd, or germline
    """
    # using a dictionary to determine 
    ibd_handler_dict: Dict = {
        "hapibd": Hapibd_Indices(ibd_program),
        "ilash": Ilash_Indices(ibd_program),
        "germine": Germline_Indices(ibd_program)
    }

    param_class = ibd_handler_dict[ibd_program.lower()]

    # updating the indices for which ever option was chosen
    param_class.update_indices()

    return param_class.return_param_dict()


def build_unique_id_dict(iid_list: list) -> Dict:
    """Function to build a dictionary of unique carriers
    Parameters
    __________
    iid_list : list
        list of iids that carry a variant based on the MEGA Ex 
        array
    
    Returns
    _______
    dict
        dictionary where the keys are the iids and the values 
        are a number
    """

    # read phenotype and build possible ID pairs
    uniqID = {}  # creates empty dictionary
    # dupID = []  # creates empty list

    IDnum = 0

    for iid in iid_list:  # This goes through each line and will get the id's

        uniqID.setdefault(iid, IDnum)

        IDnum += 1

    return uniqID

def create_ibd_arrays() -> Tuple:
    '''This creates two IBD arrays that will be used later'''

    # creating a dictionary with 22 key slots and 22 empty dictionaries
    # Also creating a dicitonary IBDindex with 22 dictionaries containing 'start': 999999999, 'end': 0, 'allpos': []
    # Using dictionary comprehension to make the two dictionaries. Just a little more concise than the for loop.
    # The 22 is for the different chromosomes.
    # the "allpos" is the breakpoints
    IBDdata = {str(i): {} for i in range(1, 23)}
    IBDindex = {
        str(i): {
            'start': 999999999,
            'end': 0,
            'allpos': []
        }
        for i in range(1, 23)
    }

    return IBDdata, IBDindex


def get_pair_string(row: pd.DataFrame, id1_indx: int, id2_indx: int, cM_indx: int, uniqID: dict) -> pd.DataFrame:
    """Function to get the pair string for each row in the dataframe
    Parameters
    __________
    row : pd.Series 
        this is the pandas series that would be each row of the chunk dataframe in the gather pairs function
    
    id1_indx : int
        this is the index to the column that has the pair1 id information
    
    id2_indx : int
        this is the index to the column that has the pair2 id information
    
    cM_indx : int
        this is the index to the column that has the total length of the ibd segment in centimorgans

    uniqID : dict 
        this is the dictionary were each key is the iids that carry the specific variant 
    
    Returns
    _______
    returns a list of pair strings
        returns the pair string 
    """
    # for row in pair_df.itertuples():

    if row[id1_indx] in (uniqID) and row[id2_indx] in (uniqID):

        if uniqID[row[id1_indx]] < uniqID[row[id2_indx]]:
                        # If both ids are in the list then it writes the pairs to a variable pair
            return '{0}:{1}-{2}'.format(row[cM_indx], row[id1_indx], row[id2_indx])

        else:
            # this just puts the ids in order
            return '{0}:{1}-{2}'.format(row[cM_indx], row[id2_indx], row[id1_indx])

    elif row[id1_indx] in (uniqID) and not row[id2_indx] in (uniqID):  # If only one id is in the uniqID then it writes it this way with the matched id in

        return '{0}:{1}-{2}'.format(row[cM_indx], row[id1_indx], row[id2_indx])

    elif row[id1_indx] not in (uniqID) and row[id2_indx] in (uniqID):  # If only id 2 is in the uniqID then it write that pair to the list
        return '{0}:{1}-{2}'.format(row[cM_indx], row[id2_indx], row[id1_indx])

def build_ibddata_and_ibddict(row: pd.Series, start_indx: int, end_indx: int, chr_indx: int, IBDdata: dict, IBDindex: dict) -> pd.Series:
    """Function that will identify breakpoints"""

    CHR: str = str(row[chr_indx])
    start: int = int(row[start_indx])
    end: int =  int(row[end_indx])
    pair: str = row["pair_string"]

    # start and end not in identified breakpoints
    if int(start) not in IBDindex[CHR]['allpos'] and int(
            end) not in IBDindex[CHR]['allpos']:

        IBDdata[CHR][str(start)] = newPOS([pair], [])
        IBDdata[CHR][str(end)] = newPOS([], [pair])
        IBDindex[CHR]['allpos'].append(int(start))
        IBDindex[CHR]['allpos'].append(int(end))

    # start is not in identified breakpoints but end is
    elif int(start) not in IBDindex[CHR]['allpos'] and int(
            end) in IBDindex[CHR]['allpos']:

        IBDdata[CHR][str(start)] = newPOS([pair], [])
        IBDdata[CHR][str(end)].rem.append(str(pair))
        IBDindex[CHR]['allpos'].append(int(start))

    # start is in identified breakpoints but end not
    elif int(start) in IBDindex[CHR]['allpos'] and int(
            end) not in IBDindex[CHR]['allpos']:

        IBDdata[CHR][str(start)].add.append(str(pair))
        IBDdata[CHR][str(end)] = newPOS([], [pair])
        IBDindex[CHR]['allpos'].append(int(end))

    # both start and end in identified breakpoints
    elif int(start) in IBDindex[CHR]['allpos'] and int(
            end) in IBDindex[CHR]['allpos']:

        IBDdata[CHR][str(start)].add.append(str(pair))
        IBDdata[CHR][str(end)].rem.append(str(pair))

    return CHR

def write_to_file(IBDdata: dict, IBDindex: dict, output: str, CHR: str, que_object, ibd_program: str, variant_name: str=None, gene_name: str=None):
    
    try:
        if len(CHR) == 1:
            chr_num: str = "0"+CHR
        else:
            chr_num = CHR

        # NEED TO FIX THIS LINE HERE
        if variant_name:
            write_path = os.path.join(
                output, "".join([ibd_program, '_', variant_name, '.chr',
                chr_num, '.small.txt.gz']))
        else: 
            write_path = os.path.join(
                output, "".join([ibd_program,'_', gene_name, '.chr',
                chr_num, '.small.txt.gz']))

        # checking to see if the file already exists from a previous 
        # one and then deleteing it
        utility_scripts.check_file(write_path)
        out = gzip.open(write_path, 'wt')

        # Writing the header line to the file
        out.write('chr\tpos\tsegments\tpairs\tadd\tdel\n')

        allibd = set([])

        for pos in sorted(IBDindex[str(CHR)]['allpos']):
            
            allibd = allibd | set(IBDdata[str(CHR)][str(pos)].add)
            allibd = allibd - set(IBDdata[str(CHR)][str(pos)].rem)

            allibdpair = {}

            if len(IBDdata[str(CHR)][str(pos)].add) == 0:
                IBDdata[str(CHR)][str(pos)].add.append('NA')
            if len(IBDdata[str(CHR)][str(pos)].rem) == 0:
                IBDdata[str(CHR)][str(pos)].rem.append('NA')

            nseg = str(len(allibd))

            for cM_pair in allibd:
                
                pair = cM_pair.split(':')[1]

                cM = cM_pair.split(':')[0]

                if pair in allibdpair:

                    allibdpair[pair] = '{0};{1}'.format(
                        allibdpair[pair], cM)

                else:

                    allibdpair[pair] = str(cM)


            npair = str(len(allibdpair))

            out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                str(CHR), str(pos), nseg, npair,
                ' '.join(IBDdata[str(CHR)][str(pos)].add),
                ' '.join(IBDdata[str(CHR)][str(pos)].rem)))

        IBDdata[str(CHR)] = []
        out.close()

        del (IBDdata)
        del (IBDindex)

    except UnboundLocalError:

        print(
            f"There were no pairs identified for the variant {variant_name}. This failure is written to a file at {''.join([output, 'nopairs_identified.txt'])}"
        )

        que_object.put(f"{variant_name}")

def fix_chr_str(chr_num: str) -> str:
    """Function that will check if the chromosome number obnly has 1 
    digit and will change the format to chrXX where X is a digit
    Parameters
    __________
    chr_num : str
        chromosome number that will either be returned with a . on 
        both sides or it will be fixed to the proper format
    
    Returns
    _______
    str
        returns a string of the correctly formatted chromosome number
    """
    chr_num = chr_num.strip(".")
    if len(chr_num) == 4:
        chr_num: str = "".join([chr_num[:3], "0", chr_num[-1]]) 
    
    return chr_num

def append_pairs(row: pd.Series, pair_info_dict: Dict, indx_dict: Dict[str, int]) -> None:
    """Function that will add the File_Pairs class to the dictionary """

    _ = pair_info_dict.setdefault("".join([row[indx_dict["id1_indx"]], "-", row[indx_dict["id2_indx"]]]), File_Pairs(row, indx_dict))

def add_pair_to_class(pair_df: pd.DataFrame, pairs_info_dict: Dict, ibd_program: str, indx_dict: Dict[str, int], chr_num: str, variant_name: str=None, gene_name: str=None) -> None:
    """Function to add class objects to a dictionary that has the IDB information for each pair
    Parameters
    __________
    pair_df : pd.DataFrame
        dataframe that is a chunk from either the iLash files or the hapibd files
    
    pairs_info_dict : Dict
        dictionary where the key will be either the variant name or the gene name and the 
        value will be another dictionary
    
    ibd_program : str
        string of either hapibd or ilash. It has the potential to use GERMLINE but that is 
        not really being used anymore
    
    indx_dict : Dict[str, int]
        Dictionary where the keys are strings that list what the index is for and the values 
        are an integer that gives the index position for that value
    
    variant_name : str
        string that list the variant name. This value is None by default
    
    gene_name : str
        string that list the gene name. This value is None by default
    """
    
    if variant_name:
        pairs_info_dict.setdefault(variant_name, {})


        pairs_info_dict[variant_name].setdefault(ibd_program, {})

        pair_df.apply(lambda row: append_pairs(row, pairs_info_dict[variant_name][ibd_program], indx_dict), axis=1)


    else:
        pairs_info_dict.setdefault(gene_name, {})

        pairs_info_dict[gene_name].setdefault(ibd_program, {})

        pair_df.apply(lambda row: append_pairs(row, pairs_info_dict[gene_name][ibd_program], indx_dict), axis=1)

    
    # need to use the apply function to add the class that has the pair info to that list
    

def gather_pairs(IBDdata: dict, IBDindex: dict, parameter_dict: dict, segment_file: str, uniqID: dict,  min_cM: int, que_object, output_path: str, ibd_program: str, pair_info_dict: Dict[str, Dict], var_position: int = None, gene_start: int = None, gene_end: int = None, variant_name=None, gene_name=None):
    '''This function will be used in the parallelism function'''
    # undoing the parameter_dict
    id1_indx = int(parameter_dict["id1_indx"])
    id2_indx = int(parameter_dict["id2_indx"])
    chr_indx = int(parameter_dict["chr_indx"])
    str_indx = int(parameter_dict["str_indx"])
    end_indx = int(parameter_dict["end_indx"])
    cM_indx = int(parameter_dict["cM_indx"])

    # getting the chromosome number that will be returned at the end of the program

    # This catches the KeyError raised because unit is only found in GERMLINE files
    try:
        unit = parameter_dict["unit"]
    except KeyError:
        unit = None

    # creating a dictionary to handle which way to filter for above the min_cM threshold
    # giving the chromosome a default value. If this value does not 
    # change throughout the program then the program will just move on
    chr_num: str = "0"

    info_dict: Dict = {}

    for chunk in pd.read_csv(segment_file,
                                sep="\t",
                                header=None,
                                chunksize=1000000):


        # Checking to see if the ids are in the uniqID dictionary
        chunk_in_uniqID: pd.DataFrame = filter_to_individual_in_uniqID(chunk, uniqID, id1_indx, id2_indx)

        
        # This is filtering the dataframe to only pairs greater than min_cM threshold
        if unit:
            # If germline files are used then you have to use this
            chunk_greater_than_3_cm: pd.DataFrame = filter_to_greater_than_3_cm(chunk_in_uniqID, cM_indx, min_cM, unit)
        else:
            chunk_greater_than_3_cm = filter_to_greater_than_3_cm(chunk_in_uniqID, cM_indx, min_cM)
        
        # filtering for values where the start value is less than the base pair and the 
        # end value is greater than the base pair
        # This method will only be done if the user is using a gene 
        # driving approach
        if var_position:
            chunk: pd.DataFrame = filter_for_correct_base_pair(chunk_greater_than_3_cm, str_indx, end_indx, var_position)

        # looking for segments where the start or end of the segment 
        # is within the gene or the 
        if gene_start and gene_end:

            chunk = filter_for_gene_site(chunk_greater_than_3_cm, str_indx, end_indx, gene_start, gene_end)
            

        if not chunk.empty:
            
            # getting rid of an warning message that indicates that the resulting chunk is a pandas dataframe
            chunk_copy: pd.DataFrame = chunk.__deepcopy__()

            # creating a new column with the pair string for each pair
            chunk_copy.loc[:,"pair_string"] = chunk_copy.apply(lambda row: get_pair_string(row, id1_indx, id2_indx, cM_indx, uniqID), axis=1)

            
            chr_num_series: pd.Series = chunk_copy.apply(lambda row: build_ibddata_and_ibddict(row, str_indx, end_indx, chr_indx, IBDdata, IBDindex), axis=1)

            chr_num: str = list(set(chr_num_series.values))[0]

            # At this point need to gather the information within a class

            # This adds the pairs string as a key to the info_dict and 
            # adds the class to that that keeps the information
            _ = chunk_copy.apply(lambda row: append_pairs(row, info_dict, parameter_dict), axis=1)
    

    if variant_name:
        _ = pair_info_dict.setdefault(variant_name, info_dict) 
           
    else:
        _ = pair_info_dict.setdefault(gene_name, info_dict)
            
    if chr_num != "0":
        if variant_name:
            write_to_file(IBDdata, IBDindex, output_path, chr_num, que_object, ibd_program, variant_name)
        if gene_name:
            write_to_file(IBDdata, IBDindex, output_path, chr_num, que_object, ibd_program, gene_name)
    else:
        if variant_name:
            print(f"the were no shared IBD segments found for the variant: {variant_name}")
        if gene_name:
            print(f"There were no shared IBD segments found for the gene {gene_name}")

    #TODO: Need to return this dictionary
            


class File_Pairs:
    """class that will keep all of the information from ilash and or hapibd"""    
    
    def __init__(self, row: pd.Series, indx_dict: Dict[str, int]) -> None:
        self.extract_information(row, indx_dict)

    def extract_information(self, row: pd.Series, indx_dict: Dict[str, int]) -> None:
        """Function to extract the information from the row series using the indices in the indx dict
        Parameter
        _________
        row : pd.Series
            pandas series that has the ibd segment information for pairs spread out over the column
            
        indx_dict : Dict[str, int]
            dictionary that has the indx names as keys and the index position as values
        """
        self.pair1: str = row[indx_dict["id1_indx"]]
        self.pair2: str = row[indx_dict["id2_indx"]]
        self.chr_num: str = row[indx_dict["chr_indx"]]
        self.str_pos: str = row[indx_dict["str_indx"]]
        self.end_pos: str = row[indx_dict["end_indx"]]
        self.cM_length: str = row[indx_dict["cM_indx"]]
        self.phase_1: str = row[indx_dict["id1_phase_indx"]]
        self.phase_2: str = row[indx_dict["id2_phase_indx"]]
