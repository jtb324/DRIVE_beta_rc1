import sys
from typing import Dict

class Parameter_Dict:
    """class to generate a dictionary of indices that certain values can be found at"""

    def __init__(self, ibd_format: str):
        self.ibd_program: str = ibd_format.lower()
        self.param_dict: dict = self.generate_base_dict()
    
    def check_program(self):
        """Function that checks if the program format is valid"""

        if self.ibd_program not in ["germline", "ilash", "hapibd"]:
            print("unrecognized or incorrect format: GERMLINE/iLASH//hap-ibd/")
            sys.exit(1)

    @staticmethod
    def generate_base_dict() -> Dict[str, int]:
        """Function to generate the base dictionary regardless of what ibd 
        program is used
        Returns
        _______
        dict
            returns a dictionary that has the indices for certain parameters 
            from the ibd files
        """
        parameter_dict = {
            "id1_indx": 0,
            "id2_indx": 2,
            "chr_indx": 4,
            "str_indx": 5,
            "end_indx": 6
        }

        return parameter_dict
    
    def return_param_dict(self) -> Dict[str, int]:
        """Function to return the parameter dictionary"""
        return self.param_dict


class Germline_Indices(Parameter_Dict):
    """class that will extend the Parameter Dict and get the extra indices"""

    def update_indices(self):
        """method to add to new indices if you use germline"""
        self.param_dict["cM_indx"] = 10
        self.param_dict["unit"] = 11

class Ilash_Indices(Parameter_Dict):
    """class that will extend the Parameter Dict and get the extra indices"""

    def update_indices(self):
        """method to add to new indices if you use ilash"""
        #this gives indexes for getting the length (the cM_indx) and the phase for each individual
        self.param_dict["cM_indx"] = 9
        self.param_dict["id1_phase_indx"] = 1
        self.param_dict["id2_phase_indx"] = 3
        

class Hapibd_Indices(Parameter_Dict):
    """class that will extend the Parameter Dict and get the extra indices"""

    def update_indices(self):
        """method to add to new indices if you use hapibd"""
        self.param_dict["cM_indx"] = 7
        self.param_dict["id1_phase_indx"] = 1
        self.param_dict["id2_phase_indx"] = 3


    