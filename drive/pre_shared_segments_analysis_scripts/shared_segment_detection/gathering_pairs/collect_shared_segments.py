#!/usr/bin/python
import sys  # THese are modules used
import gzip
import multiprocessing as mp
import re
import glob
import os
import pandas as pd
from dataclasses import dataclass
import shutil
from typing import Union

from ..generate_indx_dict.generate_dict import Germline_Indices, Ilash_Indices, Hapibd_Indices
from .filtering_functions import filter_to_greater_than_3_cm, filter_to_individual_in_uniqID, filter_for_correct_base_pair


####################################################################################################


class newPOS:
    __slots__ = 'add', 'rem'

    def __init__(self, add, rem):
        self.add = add
        self.rem = rem

class Pairs:
    """This is a class to keep information about the Pairs"""
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
    ibd_handler_dict: dict = {
        "hapibd": Hapibd_Indices(ibd_program),
        "ilash": Ilash_Indices(ibd_program),
        "germine": Germline_Indices(ibd_program)
    }

    param_class = ibd_handler_dict[ibd_program.lower()]

    # updating the indices for which ever option was chosen
    param_class.update_indices()

    return param_class.return_param_dict()


def build_unique_id_dict(iid_list: list) -> dict:
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
        # if iid in uniqID:
        #     # dupID.append(iid)
        #     pass
        # else:

        #     uniqID[iid] = IDnum
        #     IDnum = IDnum + 1

    # print('identified ' + str(len(uniqID)) + ' unique IDs')

    # Closing the file

    return uniqID

def create_ibd_arrays() -> Union[dict, dict]:
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

# class Shared_Segment_Convert(newPOS):
#     def __init__(self, shared_segment_file: str, pheno_file: str,
#                  output_path: str, ibd_program_used: str,
#                  min_cM_threshold: int, base_position,
#                  variant_id):
#         # This is the germline or hapibd or ilash file
#         self.segment_file = str(shared_segment_file)
#         # This will give the directory for where the chromosome files are found
#         self.iid_file = str(pheno_file)
#         # This will be the output directory. Need to add the ibd software to the end of it
#         self.output_dir = output_path

#         if not os.path.exists("".join([output_path, "reformatted_ibd_output/"
#                                        ])):
#             try:
#                 os.mkdir("".join([output_path, "reformatted_ibd_output/"]))
#             except FileExistsError:
#                 pass

#         self.output = "".join(
#             [output_path, "reformatted_ibd_output/", ibd_program_used])
#         self.format = str(ibd_program_used)
#         self.min_cM = int(min_cM_threshold)
#         self.bp = int(base_position)
#         # This gets the name of the variant of interest assuming it is input as a text file
#         self.variant_name = variant_id

# define a class to get some these parameters

def get_pair_string(row: pd.Series, id1_indx: int, id2_indx: int, cM_indx: int, uniqID: dict) -> pd.DataFrame:
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
    str
        returns the pair string 
    """

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

def build_ibddata_and_ibddict(row: pd.Series, start_indx: int, end_indx: int, chr_indx: str, IBDdata: dict, IBDindex: dict) -> Union[dict, dict]:
    """Function that will identify breakpoints"""
    
    CHR: str = str(row[chr_indx])
    start: int = int(row[start_indx])
    end: int =  int(row[end_indx])
    pair: str = row[len(row)-1]

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
    elif int(start_indx) in IBDindex[CHR]['allpos'] and int(
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

def gather_pairs(IBDdata: dict, IBDindex: dict, parameter_dict: dict, segment_file: str, uniqID: dict,  min_cM: int, base_pair: int):
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
        pass

    # creating a dictionary to handle which way to filter for above the min_cM threshold
    
    for chunk in pd.read_csv(segment_file,
                                sep="\t",
                                header=None,
                                chunksize=1000000):


        # Checking to see if the ids are in the uniqID dictionary
        chunk_in_uniqID = filter_to_individual_in_uniqID(chunk, uniqID, id1_indx, id2_indx)


        # This is filtering the dataframe to only pairs greater than min_cM threshold
        if unit:
            # If germline files are used then you have to use this
            chunk_greater_than_3_cm = filter_to_greater_than_3_cm(chunk_in_uniqID, cM_indx, min_cM, unit)
        else:
            chunk_greater_than_3_cm = filter_to_greater_than_3_cm(chunk_in_uniqID, cM_indx, min_cM)
        
        # filtering for values where the start value is less than the base pair and the 
        # end value is greater than the base pair
        chunk = filter_for_correct_base_pair(chunk_greater_than_3_cm, str_indx, end_indx, base_pair)

        # This will iterate through each row of the filtered chunk
        if not chunk.empty:
            chunk["pair_string"] = chunk.apply(lambda row: get_pair_string(row, id1_indx, id2_indx, cM_indx, uniqID), axis = 1)

            
            chr_num:str = chunk.apply(lambda row: build_ibddata_and_ibddict(row, str_indx, end_indx, chr_indx, IBDdata, IBDindex), axis=1)
    
    return chr_num
            

def write_to_file(IBDdata: dict, IBDindex: dict, output: str, variant_name: str, CHR: str, que_object):
    try:
        # print('identified ' + str(len(IBDindex[str(CHR)]['allpos'])) +
        #       ' breakpoints on chr' + str(CHR))

        # Opening the file .small/txt/gz file to write to
        # NEED TO FIX THIS LINE HERE
        write_path = "".join([
            output, '_', variant_name, '.chr',
            str(CHR), '.small.txt.gz'
        ])

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
    
    
