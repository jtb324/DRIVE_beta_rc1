# script for the models used during the carrier analysis
from dataclasses import dataclass
from typing import List, Optional
# from pydantic import BaseModel, validator
import models
import os
import shutil

class DirectoryNotFound(Exception):
    """Custom error message if the user inputs a directory that is not valid"""
    
    def __init__(self,message: str) -> None:
        self.message: str = message
        super().__init__(message)
@dataclass
class Carrier_Parameters():
    output: str
    pop_info_file: Optional[str] = None
    pop_code: Optional[str] = None
    run_plink: bool = True

    @classmethod
    def check_output(cls, value: str) -> None:
        """Function to make sure that the output directory exists
        Parameters
        __________
        value : str
            string that list the main directory 
            to output files to
        """
        if not os.path.isdir(value): 
            raise DirectoryNotFound(models.Colors.RED + "FATAL: "+ models.Colors.NOCOLOR + f"The provided output directory {value} was not found")
        
        return value

    def reset_output(self) -> None:
        """Function that checks if the subdirectory 'carrier_analysis_output' exists in the output directory. If it does exist then it removes it and recreates it."""

        if not os.path.isdir(os.path.join(self.output, "carrier_analysis_output")):
            print()
            os.mkdir(os.path.join(self.output, "carrier_analysis_output"))
        # If the directory exists from a previous run then it # will be removed and then an empty directory will be created
        else:
            shutil.rmtree(os.path.join(self.output, "carrier_analysis_output"))

            os.mkdir(os.path.join(self.output, "carrier_analysis_output"))

        self.full_output: str = os.path.join(self.output, "carrier_analysis_output")



