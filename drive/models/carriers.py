# script for the models used during the carrier analysis
from dataclasses import dataclass
from typing import List, Optional
import models
import os
import shutil



@dataclass
class Carrier_Analyzer:
    output: str
    pop_info_file: Optional[str] = None
    pop_code: Optional[str] = None

    def __str__(self) -> str:
        """Method to show a prettier output of the class attributes for debugging"""
        return f"output: {self.output}, population info: {self.pop_info_file}, population filter code: {self.pop_code}"

    def check_dir(self) -> None:
        """method to check if the output directory exists and if it does it will remove it and recreate the empty directoy"""

        total_directory: str = os.path.join(self.output, "carrier_analysis_output")

        if os.path.isdir(total_directory):
            shutil.rmtree(total_directory, ignore_errors=True)
        
        os.mkdir(total_directory)
        # assigns the total directory to the output
        self.output = total_directory
        
    




