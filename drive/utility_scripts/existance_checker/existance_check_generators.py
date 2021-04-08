from functools import wraps
import os

def check_dir(output_str: str, directory_name: str):
    """Function to check if the provided directory already exist and if it doesn't then it makes it
    Parameter
    _________
    output_str: str
        string listing the path output path that is being used 
        during the function that this function is run in

    directory_name : str
        the name of the directory that output will be put to 
    """
    total_directory: str = os.path.join(output_str, directory_name)

    try:
        os.mkdir(total_directory)

    except FileExistsError:
        pass

def check_file(output_str: str, file_name: str):
    """Function to check if the provided file exist from a previous run and if it does then it deletes the file
    Parameter
    _________
    output_str: str
        string listing the path output path that is being used 
        during the function that this function is run in

    file_name : str
        the name of the file
    """
    total_filepath: str = os.path.join(output_str, file_name)

    if os.path.exists(total_filepath):
        os.remove(total_filepath)

