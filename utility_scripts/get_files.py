import glob
import sys
import os

def get_file_list(files_dir: str, file_suffix: str) -> list:
    """Function to get all files in a specified directory
    Parameters
    __________
    files_dir : str
        string that describes the path to the directory that contains the
        files of interest

    file_suffix : str
        string describing the file ending that will be used to gather all 
        the different files

    Returns
    _______
    list
        returns a list containing the full file paths to all the files of 
        interest
    """
    # getting the current directory so that you can switch back to the directory that the program is running in
    cur_dir = os.getcwd()

    # switching directories to the one where the files of interest are located
    os.chdir(files_dir)

    # createing an empty list to gather the files into
    file_list = []

    
    for file in glob.glob(file_suffix):

        full_file_path = "".join([files_dir, file])

        file_list.append(full_file_path)

    os.chdir(cur_dir)

    return file_list