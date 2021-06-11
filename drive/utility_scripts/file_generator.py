# This file will contain the parent class for generating readme and log files
import os
from os import path
import getpass
import logging
from typing import Dict, Union, List
from datetime import datetime
from dataclasses import dataclass

# This is the base class for creating documentation


class Documentation:
    '''This class will be the parent class for the documentation throughout
    the program'''
    def __init__(self, file_name: str, output_path: str):

        self.file_name = "".join([output_path, file_name])
        # removing the previous file if it exist
        self.rm_previous_file()

        current_time = datetime.now()
        self.day = current_time.strftime("%b %d, %Y")
        self.time = current_time.strftime("%H:%M:%S")
        self.user = getpass.getuser()

    def rm_previous_file(self):
        '''This function will remove the readme from a previous run'''

        if path.exists(self.file_name):

            os.remove(self.file_name)

    def create_date_info(self):
        '''This function will create a file with a header'''

        with open(self.file_name, "a+") as file:

            # Writing header information to the file
            file.write(
                "################################################################\n\n"
            )
            file.write(
                "DRIVE: Distant Relatedness for Identification and Variant Evaluation\n"
            )
            file.write(f"run started on {self.day} at {self.time}\n")
            file.write(f"run started by user: {self.user}\n\n")
            file.write(
                "################################################################\n\n"
            )


class Readme(Documentation):
    '''This class will extend the Documentation class'''
    def __init__(self, file_name: str, output_path: str):

        super().__init__(file_name, output_path)

        self.write_header(output_path)

    def write_header(self, directory_name: str):
        '''This function will write a fairly genetic header for the file'''

        with open(self.file_name, "a+") as readme_file:

            readme_file.write(
                f"# README for the {directory_name} directory:\n\n")

    def add_line(self, info: Union[str, List[str]]):
        '''This function will the info_str to a new line in the document'''

        with open(self.file_name, "a+") as readme_file:

            if type(info) == str:
                readme_file.write(info + "\n")
            else:
                for line in info:
                    readme_file.write(line + "\n")

@dataclass
class Readme_Info:
    """dataclass that will hold the README body text and the readme_info
    Parameters
    __________
    readme_output_path : str
        string object that contains the filepath that the readme will be written to
     
     readme_body_text : str or List[str]
        string or list of strings that contains all the text that will be written to the readme body"""
    readme_output_path: str
    readme_body_text: Union[str, List[str]]

    readme: Readme = Readme("carrier_identification_README.md", readme_output_path)

    readme.add_line(readme_body_text)



def create_logger(output: str, name: str) -> object:
    """function to create a root logger object
    Parameters
    __________
    name : str
        name of the program being run this is the __name__ value of the initial
        program

    log_filename : str
        name of the log file that will be created

    Returns
    _______
    object
        returns a logger object that the program will log information to
    """

    # Creating a logfile with the current day as part of the log file name
    current_day = datetime.now().strftime("%m_%d_%Y")

    log_filename: str = "".join([output, current_day, "_run.log"])

    # removing a log file from a previous run
    if path.exists(log_filename):

        os.remove(log_filename)
    
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    # this next line creates a logger object at the filepath name specified
    logger = logging.getLogger(name)

    # setting the logging level to the lowest level
    logger.setLevel('DEBUG')

    # Use FileHandler() to log to a file
    file_handler = logging.FileHandler(log_filename)
    formatter = logging.Formatter(log_format)
    file_handler.setFormatter(formatter)

    # Don't forget to add the file handler
    logger.addHandler(file_handler)
    logger.info("created log file for the DRIVE analysis")
    logger.info("initializing run...")
    logger.info(f"run started by {getpass.getuser()}")

    return logger

def record_user_arguments(logger: object, user_input: Dict[str, str]):
    """Function to record the initial arguments provided by the user
    Parameters
    __________
    logger : object
        this is the logger object that information is being written to

    sys_args : object
        this is the object that contains all of the arguments passed by the user

    """
    for key, value in user_input:

        if value == None:
            logger.info(f"{key} : {value}")

def create_readmes(readme_text_list: list, output_path: str):
    """Function to create a readme
    Parameters
    __________
    readme_text_list : list
        list that contains all of the readme information

    output_path : str
        string that contains the output path to write the readme to
    """
    # creating a readme object using the Readme class
    readme = Readme("_README.md", output_path)
    # removing the previous file if it is there
    readme.rm_previous_file()
    # writing the output path to the log file
    readme.write_header(output_path)
    # creating a date tag
    readme.create_date_info()
    # adding all the text in the readme_text-list
    for readme_text in readme_text_list:
        readme.add_line(readme_text)
    
