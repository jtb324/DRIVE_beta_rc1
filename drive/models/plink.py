# models to run plink
from dataclasses import dataclass
import os
import shutil
from abc import ABC, abstractmethod
from typing import List, Optional

class RunnerInterface(ABC):

    # recode_flags: List[str] = ["recode, recodeA"]

    @abstractmethod
    def run_plink(self):
        raise NotImplementedError()

    def check_dir(self) -> None:
        """method to check if the output directory exists and if it does it will remove it and recreate the empty directoy"""

        total_directory: str = os.path.join(self.output, "plink_output_files/")

        if os.path.isdir(total_directory):
            shutil.rmtree(total_directory, ignore_errors=True)
        
        os.mkdir(total_directory)
        # assigns the total directory to the output
        self.output = total_directory

@dataclass
class Gene_Runner(RunnerInterface):
    """Class to handle running plink"""
    output: str 
    binary_file: str 
    maf_threshold: float 
    # recode_flags: List[str]
    var_file: Optional[str] = None
    

    def run_plink(self):
        """"""
        pass

@dataclass
class Range_Runner(RunnerInterface):
    """Class that runs plink if no variant file is provided"""
    output: str
    binary_file: str
    maf_threshold: float
    # recode_flags: List[str]


    def run_plink(self):
        """"""
        pass

#NOTE Need to finish these run plink methods #NOTE