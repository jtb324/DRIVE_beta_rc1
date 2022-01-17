import os
import shutil

def check_dir(output_str: str, directory_name: str) -> str:
    """Function to check if the provided directory already exist and if it doesn't then it makes it
    Parameter
    _________
    output_str: str
        string listing the path output path that is being used 
        during the function that this function is run in

    directory_name : str
        the name of the directory that output will be put to 
    
    return
    ______
    string that list the complete filepath to the directory that was made
    """
    total_directory: str = os.path.join(output_str, directory_name)

    try:
        os.mkdir(total_directory)


    except FileExistsError:
        pass

    return total_directory
     
def check_file(file_string: str):
    """Function to check if the provided file exist from a previous run and if it does then it deletes the file
    Parameter
    _________
    file_str : str
        the full path of the file
    """

    if os.path.exists(file_string):
        os.remove(file_string)

def remove_dir(directory: str) -> None:
    """Function to remove a directory that may not be 
    empty
    Parameters
    __________
    directory : str
        filepath to the directory that the user wishes 
        to remove
    """
    shutil.rmtree(directory, ignore_errors=True)