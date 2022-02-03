# This will be an object that will be used to simplify the process of determining the clusters of networks that
from dataclasses import dataclass
from abc import Abstract, abstractmethod
import os
import pandas as pd
from typing import Tuple, Dict, List
import utility_scripts

@dataclass
class Segment_Analyzer(Abstract):
    output: str

    def __post_init__(self):
        self.network_output: str = os.path.join(self.output, "networks/")
        self.segment_outputs: str = os.path.join(self.output, "formatted_ibd_output/")

    @abstractmethod
    def load_inputs(self) -> None:
        """method that will be used to load in inputs for either the phenotype search or the carrier search"""        
        pass

@dataclass
class Phenotype_Analysis(Segment_Analyzer): 

    def load_inputs(self, pheno_filepath: str, pheno_carriers: str) -> None:
        """method that loads in the two phenotype files for the phenotype analysis and loads this into 
        self attribute dictionary
        
        Parameters
        
        pheno_filepath : str
            Filepath to a text file that has the information for the gene of interest
        
        pheno_carriers : str
            Filepath to a text file that has a list of individuals that are affected by the phenotype of interest
        """

        # trying to load in the phenotype file
        pheno_df: pd.DataFrame = pd.read_csv(pheno_filepath, sep="\t", header=None)
        # trying to load in the carriers file
        carriers_df: pd.DataFrame = pd.read_csv(pheno_carriers, sep="\t")
        
        self.inputs: Dict[str, pd.DataFrame] = {
            "pheno_info": pheno_df,
            "pheno_carriers": carriers_df
            }

@dataclass
class Variant_Analysis(Segment_Analyzer):
    ibd_files: str 
    map_files_dir: str
    program_file_suffix: str
    carrier_file: str

    def load_inputs(self) -> None:
        """method that loads in files for the variant analysis and loads this into 
        self attribute dictionary
        
        Parameters
        
        
        """
        ###########################################################

        # getting the segment file list
        
        segment_file_list: List[str] = utility_scripts.get_file_list(self.ibd_files, self.program_file_suffix)
    

        # getting the list of map files
        map_file_list: List[str] = utility_scripts.get_file_list(self.map_files_dir, "*.map")

        # loading the carriers file into a dataframe
        carrier_df: pd.DataFrame = pd.read_csv(self.carrier_file, sep="\t")
        # creating a dictionary that contains the correct 
        # map, and ibd files for the correct chromosome
        # it also has a df with carriers only for each chromosome
        file_dict: dict = self.get_file_dict(map_file_list, segment_file_list, carrier_df)
        
        # filter the file dict to only the keys with values for each chromosome key
        file_dict = self.filter_empty_dictionaries(file_dict)

    def form_chr_num(self, number: int) -> str:
        """function to form the chromosome number for the key value in the get_file_dict dictionary
        Parameters
        __________
        number : int
            an integer number from the range function 
            
        Returns
        _______
        str
            returns a string of the form chrXX where X is a number
        """
        if number < 10:
            return "chr0"+str(number)
        else:
            return "chr"+str(number)

    def find_match(self, file_list: list, chr_num: str, file_type: str) -> str:
        """Function to find the file that matches the chromosome number
        Parameters
        __________
        file_list : list 
            list of files from the get_file_dict function
        chr_num : str
            chromosome number that will be used to find the match in the 
            above files   
        file_type : str
            string that indicates if it is the map file, the carrier file, or the iid file
        Returns
        _______
        str
            returns the filepath of the matched file
        """

        # if the chromosome number is less than 10 you have to remove the
        # zero for the ibd files
        if int(chr_num[-2:]) < 10:
            alt_chr_num: str = chr_num.replace("0", "")
        else: 
            alt_chr_num = chr_num
        
        # using a dictionary to do the appropriate conditioning during the 
        # list comprehension for each file type
        file_handler = {
                "map": [file for file in file_list if "".join([chr_num, "_"]) in file],
                "carrier":[file for file in file_list if "".join([chr_num, "."]) in file],
                "ibd": [file for file in file_list if "".join([alt_chr_num, "."]) in file]
                }
        
        try:
            
            file: str = file_handler[file_type][0]

        except IndexError:

            return "None"
        
        return file

    def get_file_dict(self, map_file_list: List[str], ibd_file_list: List[str], carrier_df: pd.DataFrame) -> Dict[str, Dict]:
        """Function to match up all of the files for the correct chromosome into a dictionary
        Parameters
        __________
        map_file_list : list
            a list of all the map files per chromosome as a result of running Plink

        ibd_file_list : list
            a list of all the files output by the specified ibd programs

        carrier_df : pd.Dataframe of   
            Data frame that has all the IIDs that carrier at 
            least one variant of interest. It will list the 
            variant and the chromosome for each individual
            
        Returns
        _______
        Dict[str, Dict]
            returns a dictionary where the keys are the chromsomes and the values are dictionaries with the corresponding files
        """

        # forming a dictionary of dictionaries where the key is the 
        # chromosome number
        file_dict: dict = {self.form_chr_num(i):{} for i in range(1,22)}

        
        # get the files from each list
        for key in file_dict:

            # file_dict[key].setdefault("carrier", find_match(carrier_file_list, key, "carrier"))
            file_dict[key].setdefault("map", self.find_match(map_file_list, key, "map"))
            file_dict[key].setdefault("ibd", self.find_match(ibd_file_list, key, "ibd"))

            # subset the carrier_df for just the specific chromosome 
            file_dict[key].setdefault("carrier", carrier_df[carrier_df.chr == key])


        return file_dict
    
    def filter_empty_dictionaries(self, file_dict: dict) -> dict:
        """Function to filter the file dictionary to only those keys that contain a value
        Parameters
        __________
        file_dict : dict
            dictionary contain keys for each chromosome and dictionaries 
            with filepaths as values
        
        Returns
        _______
        dict
            returns the same dictionary as inputted but only the keys with 
            values not empty dictionaries
        """

        return {key:file_dict[key] for key in file_dict if  file_dict[key]["map"] != "None"}