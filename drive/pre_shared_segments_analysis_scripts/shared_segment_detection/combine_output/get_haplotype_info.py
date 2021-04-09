import pandas as pd
from ..generate_indx_dict.generate_dict import Ilash_Indices, Hapibd_Indices

# This script will pull out the necessary information from the ibd 
# output such as start and end point, length, phase1, and phase2
class ibd_info_finder:
    """base class to get information from the ibd_files of interest"""
    def __init__(self, ibd_file: str, pair_1: str, pair_2: str, ibd_format: str) -> None:
        self.ibd_file: pd.DataFrame = pd.read_csv(ibd_file,
                                 sep="\t", header=None)
        self.pair_1: str = pair_1
        self.pair_2: str = pair_2
        self.ibd_format: str = ibd_format

class hapibd_info_finder(ibd_info_finder):
    """class to find the appropriate information for hapibd files"""
    def __init__(self, ibd_file: str, pair_1: str, pair_2: str, ibd_format: str):
        super().__init__(ibd_file, pair_1, pair_2, ibd_format)

    def get_base_indices(self):
        """Function to get the base indx dictionary from the Hapibd_Indices class"""
        # getting the base parameter indx set for each class
        hapibd_indices: Hapibd_Indices = Hapibd_Indices(self.ibd_format)
        # returning this base dict to be a class parameter
        self.indx_dict: dict = hapibd_indices.return_param_dict()
    
    def update_base_indx_dict(self):
        """Function to update the indx_dict so that it has the appropriate indices"""
        # generating the index for the cM total length
        self.indx_dict["cM_indx"] = 7
        # getting the phasing for pair 1
        self.indx_dict["phase_1"] = 1
        # getting the phasing for pair 2
        self.indx_dict["phase_2"] = 3


class ilash_info_finder(ibd_info_finder):
    """class to find the appropriate information for ilash files"""
    def __init__(self, ibd_file: str, pair_1: str, pair_2: str, ibd_format: str):
        super().__init__(ibd_file, pair_1, pair_2, ibd_format)

    def get_base_indices(self):
        """Function to get the base indx dictionary from the Ilash_Indices class"""
        # getting the base parameter indx set for each class
        ilash_indices: Ilash_Indices = Ilash_Indices(self.ibd_format)
        # returning this base dict to be a class parameter
        self.indx_dict: dict = ilash_indices.return_param_dict()

    def update_base_indx_dict(self):
        """Function to update the indx_dict so that it has the appropriate indices"""
        # generating the index for the cM total length
        self.indx_dict["cM_indx"] = 9
        # getting the phasing for pair 1
        self.indx_dict["phase_1"] = 1
        # getting the phasing for pair 2
        self.indx_dict["phase_2"] = 3
