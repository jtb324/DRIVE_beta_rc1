from pydantic import BaseModel, validator, root_validator
from typing import List, Dict, Optional
import os
import models

class IncorrectIBDProgramError(Exception):
    """Class that handles the error if the user supplies a unsupported ibd program"""
    def __init__(self, ibd_list: List[str], message: str) -> None:
        self.ibd_programs: List[str] = ibd_list
        self.message: str = message
        super().__init__(message)

class IncorrectMafThresholdError(Exception):
    """Custom class that handles the error if the user provides an incorrect maf threshold"""
    def __init__(self, threshold: float, message: str) -> None:
        self.threshold: float = threshold
        self.message: str = message
        super().__init__(message)

class DirectoryNotFoundError(Exception):
    """Custom error message if the user inputs a directory that is not valid"""
    
    def __init__(self, directory: str, message: str) -> None:
        self.directory: str = directory
        self.message: str = message
        super().__init__(message)

class IllogicalGeneRangeError(Exception):
    """Custom error message if the gene range is illogical"""
    def __init__(self, gene_range: List[str], message: str) -> None:
        self.gene_range: List[str] = gene_range
        self.message: str = message
        super().__init__(message)

class ImproperPopulationCode(Exception):
    """Custom error message if the population code is not in an accepted format"""
    def __init__(self, population_code: str, message: str) -> None:
        self.population_code: str = population_code
        self.message: str = message
        super().__init__(message)

class UnsupportedFileType(Exception):
    """Custom error message if the variant file has any file extension except csv, xlsx, or '' """
    def __init__(self, extension: str, message: str) -> None:
        self.extension: str =extension
        self.message: str = message
        super().__init__(message)

class InputParams(BaseModel):
    title: str 
    output_path: str
    ethnicity_path: Optional[str]
    pop_code: str 
    pop_info_file: str
    plink_parameters: Dict
    var_file: str
    bfile: str
    pfiles: Optional[str]
    carrier_file_dir: str
    pheno_gmap: str 
    pheno_carriers: str
    gene_range: List[int]
    ibd_programs: List[str]
    min_cm: int
    maf_threshold: float

    def __init__(self, **kwargs) -> None:
        """Class that takes in the user configurations for the program"""
        # These next sections map out the specific inner files of the toml file into the 
        # corresponding model attributes
        kwargs["output_path"] = kwargs["output"]["output"]

        kwargs["ethinicity_path"] = kwargs["population_parameters"]["ethnicity_file"]

        kwargs["pop_code"] = kwargs["population_parameters"]["pop_code"]

        kwargs["pop_info_file"] = kwargs["population_parameters"]["pop_info_file"]

        kwargs["var_file"] = kwargs["inputs"]["variant_file"]

        kwargs["bfile"] = kwargs["inputs"]["binary_file"]

        kwargs["pfile"] = kwargs["inputs"]["ped_files"]

        kwargs["carrier_file_dir"] = kwargs["inputs"]["carrier_files"]

        kwargs["pheno_gmap"] = kwargs["inputs"]["pheno_gene_info"]

        kwargs["pheno_carriers"] = kwargs["inputs"]["phenotype_carriers"]

        kwargs["min_cm"] = kwargs["thresholds"]["min_cm"]

        kwargs["maf_threshold"] = kwargs["thresholds"]["maf_threshold"]

        kwargs["gene_range"] = kwargs["gene_range"]["gene_range"]

        kwargs["ibd_programs"] = kwargs["ibd_programs"]["ibd_programs"]

        super().__init__(**kwargs)

    @validator('ibd_programs')
    @classmethod
    def incorrect_ibd_program(cls, value: List[str]) -> List[str]:
        """Function to check if the user has provided an incorrect ibd_program"""

        # The program supports ilash, hapibd, and germline
        supported_ibd_programs: List[str] = ["hapibd", "ilash", "germline"]

        # checking to see if the supplied programs are in the above list. If they are
        # not then an error is raised that tells the user what the acceptable 
        # programs are
        if len([program for program in value if program.lower() not in supported_ibd_programs]) != 0:
            error_message: str = models.Colors.RED.value + "FATAL: " + models.Colors.NOCOLOR.value + f"One or more of the ibd programs you supplied are not supported. Supported programs are {', '.join(supported_ibd_programs)}\n"
            raise IncorrectIBDProgramError(value, error_message)
        return value

    @validator('maf_threshold')
    @classmethod
    def incorrect_maf_threshold(cls, value: float) -> float:
        """Function that makes sure that the maf_threshold value is >0 and <=1"""

        # if the threshold is incorrect than it raises an error
        if value <= 0 or value >= 0.5:
            error_message: str = models.Colors.RED.value + "FATAL: " + models.Colors.NOCOLOR.value +f"Incorrect minor allele frequency threshold supplied. The minor allele frequency threshold is expected to be > 0 or < 0.5"
            raise IncorrectMafThresholdError(value, error_message)
        return value
    
    @validator('output_path')
    @classmethod
    def output_path_not_found(cls, value: str) -> str:
        """Function to check and make sure that the output path exists"""
        # checks if the specified directory exists. If it doesn't then it raises 
        # an error_message
        if not os.path.exists(value):
            error_message: str = models.Colors.RED.value + "FATAL: "+ models.Colors.NOCOLOR.value + f"The provided output directory {value} was not found"
            raise DirectoryNotFoundError(error_message)
        return value
    
    @validator('gene_range')
    @classmethod
    def illogical_gene_range(cls, value: List[int]) -> List[str]:
        """Function to make sure that the 1st value in the gene range is smaller than the second value"""

        # Checks if the first gene range value is larger than the second value which 
        # would indicate that the gene range is reversed
        if value[0] > value[1]:
            error_message: str = models.Colors.RED.value + "FATAL: "+ models.Colors.NOCOLOR.value + f"The gene range provided was illogical. The second value was smaller than the first value"
            raise IllogicalGeneRangeError(value, error_message)

    @validator('pop_code')
    @classmethod
    def nonsupported_pop_code(cls, value: str) -> str:
        """Function to check if the population code is supported"""

        pop_codes: List[str] = ["", "EUR", "EAS", "AFR", "AMR", "SAS"]

        if value not in pop_codes:
            error_message: str = models.Colors.RED.value + "FATAL: "+ models.Colors.NOCOLOR.value + f"The population code provided was not a valid value, excepted values are: ''. {', '.join(pop_codes)}"
            raise ImproperPopulationCode(value, error_message)

    @validator("var_file")
    @classmethod
    def incorrect_variant_file(cls, value: str) -> str:
        """Function to check if the extension of the var file is supported"""
        extension_types: List[str] = ["", ".csv", ".xlsx"]

        # determining where the last '.' is in the filename 
        # that way we can get the extension. If there is no 
        # '.' then an IndexError is raised and the value is 
        # an empty string. 
        try:
            # find the index position of the last dot
            last_dot_indx: int = [i for i, letter in enumerate(value) if letter == "."][-1]

            extension: str = value[last_dot_indx:]

        except IndexError:

            extension: str = ""

        print(f"extension: {extension}")
        if extension not in extension_types:
            error_message: str = models.Colors.RED.value + "FATAL: "+ models.Colors.NOCOLOR.value + f"The extension of the variant files is not correct. If the file is not provided then the configuration file should be ''. If the file is provided than the extension needs to be either .xlsx or .csv"
            raise UnsupportedFileType(extension, error_message)

