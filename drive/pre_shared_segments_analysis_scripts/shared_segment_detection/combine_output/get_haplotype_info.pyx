import pandas as pd
import gzip
from typing import List, Dict, Tuple
from ..generate_indx_dict.generate_dict import Ilash_Indices, Hapibd_Indices

# This script will pull out the necessary information from the ibd 
# output such as start and end point, length, phase1, and phase2
class ibd_info_finder:
    """base class to get information from the ibd_files of interest"""
    def __init__(self,ibd_file: pd.DataFrame,str pair_1,str pair_2, str ibd_format) -> None:
        """base class that will be used to get information about the ibd segment.
        Parameters
        __________
        ibd_df : pd.DataFrame
            dataframe of the output from either ilash or hapibd
        
        pair_1 : str
            string of the first iid in the pair
        
        pair_2 : str
            string of the second iid in the pair
        
        ibd_format : str
            string of the ibd program that was used in the analysis. This will be 
            either ilash, hapibd, or both of the programs
        """
        self.ibd_file: str = ibd_file
        self.pair_1: str = pair_1
        self.pair_2: str = pair_2
        self.ibd_format: str = ibd_format
        self.indx_dict: Dict[str, int] = {}
    
    def get_len_info(self, var_position: int=None, gene_start: int =None, gene_end: int = None) -> dict:
        """Function to get the shared segment lengths for each pair
        Parameters
        __________
        var_position : int
            This is the base position of the variant of interest. This value is None 
            by default
        
        gene_start : int
            This is the base position of where the gene of interest starts. This 
            value is None by default.
            
        gene_end : int
            This is the base position of where the gene of interest ends. This value is None by default.
        
        Returns
        _______
        dict
            returns a dictionary that has the start, end, length, phase1, and phase2 
            information for the pair of interest
        """

        # filtering the ibd file for just values that have the pair 1 or pair 2
        filtered_df: pd.DataFrame = self.ibd_file[(self.ibd_file[self.indx_dict["id1_indx"]] == self.pair_1) & (self.ibd_file[self.indx_dict["id2_indx"]] == self.pair_2)]

        # now need to filter to just the values that have are have 
        # the var_position within them if the phenotype analysis is 
        # not choosen or the values that involve the gene start and 
        # end
        # doing the non phenotype analysis first
        # makes sure the variant postion is within the gene

        if filtered_df.empty:
            
            # switching the pairs if the dataframe is empty
            filtered_df: pd.DataFrame = self.ibd_file[(self.ibd_file[self.indx_dict["id1_indx"]] == self.pair_2) & (self.ibd_file[self.indx_dict["id2_indx"]] == self.pair_1)]

            # if it is still empty than it returns a null dictionary
            if filtered_df.empty:
                return {
                    "start": "N/A",
                    "end": "N/A",
                    "length": "N/A",
                    "phase1": "N/A",
                    "phase2": "N/A"
                }
        # if the var_position is provided then the function will 
        # filter for segments that have the variant within them
        if var_position:
            within_segment_df: pd.DataFrame = filtered_df[(filtered_df[self.indx_dict["str_indx"]] <= var_position) & (filtered_df[self.indx_dict["end_indx"]] >= var_position)]
        # otherwise it will filtered based on the gene start and end 
        # position. This is only uised for the phenotype analysis
        else:
            within_segment_df: pd.DataFrame = filtered_df[((filtered_df[self.indx_dict["end_indx"]].astype(int) <= int(gene_end)) & (filtered_df[self.indx_dict["end_indx"]].astype(int) >= int(gene_start)))| ((filtered_df[self.indx_dict["str_indx"]].astype(int) >= int(gene_start)) & (filtered_df[self.indx_dict["str_indx"]].astype(int) <= int(gene_end))) | ((filtered_df[self.indx_dict["str_indx"]].astype(int) <= int(gene_start)) & (filtered_df[self.indx_dict["end_indx"]].astype(int) >= int(gene_end)))]
        
        # if the dataframe is empty then return a dictionary of just 
        # N/A's
        if within_segment_df.empty:
            return {
                "start": "N/A",
                "end": "N/A",
                "length": "N/A",
                "phase1": "N/A",
                "phase2": "N/A"
            }
        else:
            # assigning the phase information to a variable
            phase_1: str = within_segment_df[self.indx_dict["phase_1"]].values.tolist()[0]
            phase_2: str = within_segment_df[self.indx_dict["phase_2"]].values.tolist()[0]

            # For ilash the phase is either 0 or 1 and is at the end of the string 
            # so this line checks if the phase string length is greater than 1 and 
            # if it is it will get the 0 or one from the end of the string. This 
            # step is not necessary for hapibd because the phase is just a single 
            # digit of 1 or 2
            if len(str(phase_1)) != 1:
                phase_1 = phase_1[-1]
                phase_2 = phase_2[-1]

            return {
                "start": within_segment_df[self.indx_dict["str_indx"]].values.tolist()[0],
                "end": within_segment_df[self.indx_dict["end_indx"]].values.tolist()[0],
                "length": within_segment_df[self.indx_dict["cM_indx"]].values.tolist()[0],
                "phase1": phase_1,
                "phase2": phase_2
            }

class hapibd_info_finder(ibd_info_finder):
    """class to find the appropriate information for hapibd files"""
    def __init__(self, ibd_file: pd.DataFrame, str pair_1: str, str pair_2: str, str ibd_format: str):
        super().__init__(ibd_file, pair_1, pair_2, ibd_format)
        self.get_base_indices()

    def get_base_indices(self):
        """Function to get the base indx dictionary from the Hapibd_Indices class. 
        Assigns this base indx dictionary to the attribute indx_dict"""
        # getting the base parameter indx set for each class
        hapibd_indices: Hapibd_Indices = Hapibd_Indices(self.ibd_format)
        # returning this base dict to be a class parameter
        self.indx_dict: dict = hapibd_indices.return_param_dict()

        # generating the index for the cM total length
        self.indx_dict["cM_indx"] = 7
        # getting the phasing for pair 1
        self.indx_dict["phase_1"] = 1
        # getting the phasing for pair 2
        self.indx_dict["phase_2"] = 3
    
        


class ilash_info_finder(ibd_info_finder):
    """class to find the appropriate information for ilash files"""
    def __init__(self, ibd_file: pd.DataFrame, str pair_1, str pair_2, str ibd_format):
        super().__init__(ibd_file, pair_1, pair_2, ibd_format)
        self.get_base_indices()

    def get_base_indices(self):
        """Function to get the base indx dictionary from the Ilash_Indices class
        Assigns this base indx dictionary to the attribute indx_dict"""
        # getting the base parameter indx set for each class
        ilash_indices: Ilash_Indices = Ilash_Indices(self.ibd_format)
        # returning this base dict to be a class parameter
        self.indx_dict: dict = ilash_indices.return_param_dict()
        # generating the index for the cM total length
        self.indx_dict["cM_indx"] = 9
        # getting the phasing for pair 1
        self.indx_dict["phase_1"] = 1
        # getting the phasing for pair 2
        self.indx_dict["phase_2"] = 3

        
