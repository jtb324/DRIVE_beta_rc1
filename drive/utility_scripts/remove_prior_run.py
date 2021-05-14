import os
import shutil
from typing import List
import glob
from functools import wraps

# This decorator will delete directories and log files from the previous run and will then recreate the necessary directories
def check_file_existance(func):
    @wraps(func)
    def inner_decorator(output: str):
        """Decorator that will be used to check if the output from a prior run is there, It will delete these files if they are there 
        Parameters
        __________
        output : str
            string that will be the directory that files will be put out into
        """
        file_directory_list: List[str] = [
            "carrier_analysis_output",
            "formatted_ibd_output",
            "plink_output_files",
        ]
        # have to form each directory path
        for file in file_directory_list:
            directory_path: str = os.path.join(output, file)

            func(directory_path)

            make_directory(directory_path)

        # removing all the log files from a previous run
        remove_log_files(output)

        
    return inner_decorator

@check_file_existance
def remove_dir(output_dir: str):
    """function to remove all the output directories from any previous run
    Parameters
    __________
    output_dir : str 
        string that list the output directory that all the output folders will be written into 
    """
    shutil.rmtree(output_dir, ignore_errors=True)

def remove_log_files(output_dir: str):
    """function that will remove the log files
    Parameters
    __________
    output_dir : str 
        string that list the output directory that all the output folders will be written into
    """
    current_dir: str = os.getcwd()

    os.chdir(output_dir)

    for file in glob.glob("*.log"):

        os.remove(os.path.join(output_dir, file))

    os.chdir(current_dir)

def make_directory(output_dir: str):
    """function that will create necessary output directories
    Parameters
    __________
    output_dir : str 
        string that list the output directory that all the output folders will be written into
    """
    os.mkdir(output_dir)