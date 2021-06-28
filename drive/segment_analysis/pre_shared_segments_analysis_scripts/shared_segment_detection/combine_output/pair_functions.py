import pandas as pd
from .get_haplotype_info import hapibd_info_finder, ilash_info_finder

import sys
from typing import List, Dict, Optional
# This script keeps some of the functions that are used for determining if pairs are found

def is_max_pairs_found(curr_max_pairs: int, new_max_pairs: int) -> int:
    """Function to determine if the max number of pairs was found.
    It checks to see if the max pair from the previous row is largeer or not
    Parameters
    __________
    curr_max_pairs : int
        integer of the number of max pairs from the previous row
        
    new_max_pairs : int
        integer of the number of max pairs found for the current row
        
    Returns
    _______
    int
        integer of either 1 or 0 where 1 means that the max pair was found in the previous row and 0 if the new max pairs is larger than the previous row
    """

    pair_handler_dict = {False: 0, True: 1}

    # This function will return either 1 or zero from the pair handler_dict based on whether or not curr_max_pair from the previous row is less than or greater than the max pairs from the current row
    return pair_handler_dict[curr_max_pairs >= new_max_pairs]

def after_max_pair_found(curr_max_pair: int, new_max_pair: int) -> int:
    '''This function will break out of the loop if a max pair is found and then the size of the pair list starts decreasing'''

    # If the previous rows max number of pairs is higher than the current row then this function will return a one
    if curr_max_pair > new_max_pair:
        return 1
    elif curr_max_pair == new_max_pair:  # If the the above is not true than it returns a 0
        return 0

# creating a class to keep information about the pair string
class Pairs:
    """Class to get specific information about the pairs such as 
    The ibd_programs that identified the pairs, the pair1, pair2, and then making a
    
    Parameters
    __________
    pair_str : str
        string that has all the pairs with there ibd program
    """
    def __init__(self, pair_str: str):
        self.program: str = pair_str.split(":")[0]
        self.pair1: str = pair_str.split(":")[1].split("-")[0].strip("\n")
        self.pair2: str = pair_str.split(":")[1].split("-")[1].strip("\n")


# class that will keep track of certain information for the max Pair row
class Pair_Info_Class:
    """class to organize information about the pair 1 and pair 2 
    Parameters
    __________
    pair_string : str
        string containing the grid ID for Pair 1 and pair 2

    identifier : str
        the variant that the pairs are being tested for
    
    chromo_num : str
        string describing the chromosome number that the shared 
        segment is on
    allpair_filepath : str
        string containing the filepath that the allpair file will be output to
    """

    def __init__(self, pair_string: str, identifier: str, chromo_num: str, allpair_filepath: str, analysis_type: str):
        self.pair_list: list = self.split_pair_str(pair_string)
        self.identifier: str = identifier
        self.chromo_num: str = chromo_num
        self.output_path: str = allpair_filepath
        self.analysis_type: str = analysis_type

    
    @staticmethod
    def split_pair_str(pair_string: str) -> List[str]:
        """Function to split the string of pairs 
        Parameters
        __________
        pair_string : str
            string containing the grid ids for pair 1 and pair 2
        Returns
        _______
        list 
            list containing the id for pair 1 and pair 2
        """
        # getting the string of just the two pairs
        pairs_str: str = pair_string.split("\t")[4]

        # splitting the above string into a list with pair 1 and
        # pair 2
        pairs_list: list = pairs_str.split(" ")

        return pairs_list

    def iid_list_handler(self, carrier_df: str = None, pheno_carriers: pd.DataFrame = None):
        """Function to handle which iid list method will be formed.

        Parameters
        __________
        carrier_dir : str
            directory that specifies where the 
            single_variant_carrier.csv files are that list which
            iid carries which variant

        pheno_carriers : pd.DataFrame
            dataframe from the pheno_carriers.xlsx file that has
            two columns where one is the iid and the other is the 
            gene. This file contains iids who were identified as
            having the phenotype and then the gene that the 
            phenotype is associated with in the next column
        """
        
        if self.analysis_type == "phenotype":
            self.carrier_iid_list(pheno_carriers=pheno_carriers)
        else:
            self.carrier_iid_list(carrier_analysis_df =carrier_df)
        # iid_option_handler: dict = {
        #     True: self.carrier_iid_list(pheno_carriers=pheno_carriers), 
        #     False: self.carrier_iid_list(carrier_analysis_dir=carrier_dir)
        # }

        # iid_option_handler[self.analysis_type == "phenotype"]

    def carrier_iid_list(self, carrier_analysis_df: Optional[pd.DataFrame] = None, pheno_carriers: Optional[pd.DataFrame] = None):
        """Function to determine the iid list based on the analysis 
        method. Assigns this list to an attribute called iid_list of 
        the Pair_Info_Class object   

        Parameters
        __________
        carrier_analysis_dir : Optional[pd_dataframe]
            directory that specifies where the 
            single_variant_carrier.csv files are that list which
            iid carries which variant

        pheno_carriers : pd.DataFrame
            dataframe from the pheno_carriers.xlsx file that has
            two columns where one is the iid and the other is the 
            gene. This file contains iids who were identified as
            having the phenotype and then the gene that the 
            phenotype is associated with in the next column 
        """
        
        # if it is a phenotype approach then we need to get all the 
        # carriers for a specific gene

        if self.analysis_type == "phenotype":
            
            # pulling out the iid_list from the pheno carriers for the specific gene that you have
            iid_list: List[str] = pheno_carriers[pheno_carriers["gene"] == self.identifier]["IID"].values.tolist()


        # If the analysis is not the phenotype analysis then need to 
        # gather the carrier files from the specified directory
        # This will find individuals that carry a specific variant
        else:

            # filter dataframe for those individuals carrying the variant
            iid_list: List[str] = carrier_analysis_df["iid"].values.tolist()

        self.iid_list: List[str] = iid_list
    
    def check_for_missed_carriers(self, pair_object: Pairs) -> int:
        """This function will check for variants that may be missed carriers, meaning that pair 2 may be connected to other carriers

        Parameters
        __________
        pair_object : Pairs
            Pairs class object that has information about the pair1 and pair2 grid 
            iid and then the ibd programs that identified the IBD segment for the 
            pair
        
        Returns
        _______
        int
            returns an integer that is the number of carriers that the pair 2 is connected to 
        """
        # get a list of all strings that contain the pair2 value
        pair2_in_str: list = [string for string in self.pair_list if pair_object.pair2 in string]

        # getting a list of connected carriers that can be summed at the end
        connected_carriers_list: list = []

        # iterating through each pair
        for pair in pair2_in_str:

            pair_iids: list = pair.split(":")[1]

            pair_1, pair_2 = pair_iids.split("-")

            # checking if the pair_1 and pair_2 are a carrier
            if pair_1 in self.iid_list or pair_2 in self.iid_list:

                connected_carriers_list.append(pair)

        return int(len(connected_carriers_list))

    def set_missed_carrier_status_null(self, pairs_dict: dict, pairs_object: Pairs):
        """Function to provide values to the dictionary when the pair2 is a confirmed carrier and then the connected carriers value doesn't matter
        
        Parameters
        __________
        pairs_dict : dict
            dictionary that contains information about the pairs that are being analyzed
        
        pairs_object : Pairs
            pairs class object that contains information about the pairs being 
            analyzed
        """
        # setting the connected_carriers value to N/A because these values can't be 
        # missed carriers
        pairs_dict[(pairs_object.pair1, pairs_object.pair2)]["connected_carriers"] = "N/A"
        # setting the missed carrier value to 0 because the pair can't be a missed 
        # carrier since it is a confirmed carrier
        pairs_dict[(pairs_object.pair1, pairs_object.pair2)]["missed_carrier"] = 0

    def set_missed_carrier_status(self, pairs_dict: dict, connected_carriers: int, pairs_object: Pairs):
        """Function to determine if the pair_2 may be a missed carrier based on the 
        number of connected_carriers
        
        Parameters
        __________
        pairs_dict : dict
            dictionary that contains information about the pairs that are being analyzed
        
        connected_carriers : int
            integer of the number of connected carriers that pair two is connected to
        
        pairs_object : Pairs
            pairs class object that contains information about the pairs being 
            analyzed
        """
        # setting the connected_carriers values 
        pairs_dict[(pairs_object.pair1, pairs_object.pair2)]["connected_carriers"] = connected_carriers

        if connected_carriers >= 2:
            pairs_dict[(pairs_object.pair1, pairs_object.pair2)]["missed_carrier"] = 1
        else:
            pairs_dict[(pairs_object.pair1, pairs_object.pair2)]["missed_carrier"] = 0
        
    def generate_pairs_dict(self, analysis_type_info: Dict, pair_info_dict: Dict[str, Dict], identifier) -> List[str]:
        """Function that generates a dictionary of all the pairs and information
        
        Parameters
        __________
        hapibd_file : pd.DataFrame
            dataframe of the hapibd output file for the specific chromosome
        
        ilash_file : pd.DataFrame
            dataframe of the ilash output file for the specific chromosome
        
        analysis_type_info : dict
            dictionary containing the analysis type, the variant position or the gene start and end. These will be under these respective keys

        pair_info_dict : Dict[str, Dict[str, Dict[str, object]]]
            dictionary of dictionaries where the outer keys are either 
            ilash or hapibd and the value is a dictionary where the key 
            is the MEGA probe id, and the value is a dictionary where 
            the key is the pair in the form of 'pair1-pair2' and the 
            value is a class object that has information

        Returns
        _______
        list 
            returns a list that has information about the for each 
            pair as the values"""
        # creating an empty dictionary to put the pairs into
        pairs_dict: dict = {}

        pairs_list: list = []
        # iterating through each pair in the pair list
        for pair in self.pair_list:
            
            
            # creating an Pairs object that has each pair string
            pair_object: Pairs = Pairs(pair)

            # creating a string of the pairs to identify the pairs from 
            # the dictionary
            pair_str: str = "".join([pair_object.pair1,"-", pair_object.pair2])

            # creating an alternate string to check if the pairs are in the reverse order
            alt_pair_str: str = "".join([pair_object.pair2,"-", pair_object.pair1])

            # updating the pairs_dict so that there is a key for this pair
            pairs_dict[(pair_object.pair1, pair_object.pair2)] = {}

            # check if the second pair is in the carrier iid list
    
            
            if pair_object.pair2 in self.iid_list:

                
                # if the pair is a carrier than the carrier status will be set to 1
                pairs_dict[(pair_object.pair1, pair_object.pair2)]["carrier_status"] = 1

                # assigning values for the missed carriers
                self.set_missed_carrier_status_null(pairs_dict, pair_object)

                 # need to 
            else:
                # if the pair is not carrier than the carrier status will be set to 0
                pairs_dict[(pair_object.pair1, pair_object.pair2)]["carrier_status"] = 0

                # determining if that pair may be connected to other carriers
                connected_carriers: int = self.check_for_missed_carriers(pair_object)

                # assigning values for the missed carriers
                self.set_missed_carrier_status(pairs_dict, connected_carriers, pair_object)

            # # need to get the information from hapibd and ilash about the segment length
            # hapibd_info_object: hapibd_info_finder = hapibd_info_finder(hapibd_file, pair_object.pair1, pair_object.pair2, "hapibd")

            # ilash_info_object: ilash_info_finder = ilash_info_finder(ilash_file, pair_object.pair1, pair_object.pair2, "ilash")
            
            # getting the ilash and hapibd objects with all the information for the pair

            try:
                hapibd_pair_object = pair_info_dict["hapibd"][identifier][pair_str]
            
            except KeyError:
                hapibd_pair_object = pair_info_dict["hapibd"][identifier].get(alt_pair_str)

            try:
                ilash_pair_object = pair_info_dict["ilash"][identifier][pair_str]

            except KeyError:
                ilash_pair_object = pair_info_dict["ilash"][identifier].get(alt_pair_str)

            
            # need to create an object that has can differiate the objects based on if hapibd_pair_object is None or if ilash_pair_object is none
            if hapibd_pair_object and ilash_pair_object:
            # if the analysis type is phenotype then have to use the 
            # gene start and end
                if analysis_type_info["analysis_type"] == "phenotype":
                    
                    # hapibd_info_dict: dict = hapibd_info_object.get_len_info(gene_start=analysis_type_info["gene_start"], gene_end=analysis_type_info["gene_end"])

                    # ilash_info_dict: dict = ilash_info_object.get_len_info(gene_start=analysis_type_info["gene_start"], gene_end=analysis_type_info["gene_end"])

                    # creating the pair_str when the phenotype analysis is selected
                    pair_str:str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{'N/A'}\t{self.identifier}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\t{hapibd_pair_object.phase_1}\t{hapibd_pair_object.phase_2}\t{ilash_pair_object.phase_1}\t{ilash_pair_object.phase_2}\t{hapibd_pair_object.str_pos}\t{hapibd_pair_object.end_pos}\t{hapibd_pair_object.cM_length}\t{ilash_pair_object.str_pos}\t{ilash_pair_object.end_pos}\t{ilash_pair_object.cM_length}\n"

                    # pair_str:str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{'N/A'}\t{self.identifier}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\n"
                # if it is not phenotype then you have to use the 
                # variant position
                else:
                    # hapibd_info_dict: dict = hapibd_info_object.get_len_info(var_position=analysis_type_info["variant_pos"])

                    # ilash_info_dict: dict = ilash_info_object.get_len_info(var_position=analysis_type_info["variant_pos"])
                
                    # writing the pair string when the phenotype analysis is not used
                    pair_str: str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{self.identifier}\t{'N/A'}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\t{hapibd_pair_object.phase_1}\t{hapibd_pair_object.phase_2}\t{ilash_pair_object.phase_1}\t{ilash_pair_object.phase_2}\t{hapibd_pair_object.str_pos}\t{hapibd_pair_object.end_pos}\t{hapibd_pair_object.cM_length}\t{ilash_pair_object.str_pos}\t{ilash_pair_object.end_pos}\t{ilash_pair_object.cM_length}\n"

                    # pair_str: str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{self.identifier}\t{'N/A'}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\n"
            elif not hapibd_pair_object and not ilash_pair_object:

                if analysis_type_info["analysis_type"] == "phenotype":

                    # creating the pair_str when the phenotype analysis is selected
                    pair_str:str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{'N/A'}\t{self.identifier}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\n"

                # if it is not phenotype then you have to use the 
                # variant position
                else:
                    
                    # writing the pair string when the phenotype analysis is not used
                    pair_str: str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{self.identifier}\t{'N/A'}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{'N/A'}\n"

            elif not hapibd_info_finder:
                
            # if the analysis type is phenotype then have to use the 
            # gene start and end
                if analysis_type_info["analysis_type"] == "phenotype":
                    
                    
                    # creating the pair_str when the phenotype analysis is selected
                    pair_str:str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{'N/A'}\t{self.identifier}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\t{'N/A'}\t{'N/A'}\t{ilash_pair_object.phase_1}\t{ilash_pair_object.phase_2}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{ilash_pair_object.str_pos}\t{ilash_pair_object.end_pos}\t{ilash_pair_object.cM_length}\n"

                    
                # if it is not phenotype then you have to use the 
                # variant position
                else:
                    
                
                    # writing the pair string when the phenotype analysis is not used
                    pair_str: str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{self.identifier}\t{'N/A'}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\t{'N/A'}\t{'N/A'}\t{ilash_pair_object.phase_1}\t{ilash_pair_object.phase_2}\t{'N/A'}\t{'N/A'}\t{'N/A'}\t{ilash_pair_object.str_pos}\t{ilash_pair_object.end_pos}\t{ilash_pair_object.cM_length}\n"

            elif not ilash_pair_object:

            # if the analysis type is phenotype then have to use the 
            # gene start and end
                if analysis_type_info["analysis_type"] == "phenotype":

                    # creating the pair_str when the phenotype analysis is selected
                    pair_str:str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{'N/A'}\t{self.identifier}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\t{hapibd_pair_object.phase_1}\t{hapibd_pair_object.phase_2}\t{'N/A'}\t{'N/A'}\t{hapibd_pair_object.str_pos}\t{hapibd_pair_object.end_pos}\t{hapibd_pair_object.cM_length}\t{'N/A'}\t{'N/A'}\t{'N/A'}\n"

                # if it is not phenotype then you have to use the 
                # variant position
                else:
                    
                    # writing the pair string when the phenotype analysis is not used
                    pair_str: str = f"{pair_object.program}\t{pair_object.pair1}\t{pair_object.pair2}\t{self.chromo_num.strip('.')[3:]}\t{self.identifier}\t{'N/A'}\t{pairs_dict[(pair_object.pair1, pair_object.pair2)]['carrier_status']}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['missed_carrier'])}\t{str(pairs_dict[(pair_object.pair1, pair_object.pair2)]['connected_carriers'])}\t{hapibd_pair_object.phase_1}\t{hapibd_pair_object.phase_2}\t{'N/A'}\t{'N/A'}\t{hapibd_pair_object.str_pos}\t{hapibd_pair_object.end_pos}\t{hapibd_pair_object.cM_length}\t{'N/A'}\t{'N/A'}\t{'N/A'}\n"

            

            pairs_list.append(pair_str)

        return pairs_list


