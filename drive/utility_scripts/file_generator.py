# This file will contain the parent class for generating readme and log files
import os
from os import path
import getpass
import logging
from typing import Dict
from datetime import datetime

# This is the base class for creating documentation


class Documentation:
    '''This class will be the parent class for the documentation throughout
    the program'''
    def __init__(self, file_name: str, output_path: str):

        self.file_name = "".join([output_path, file_name])
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

    def write_header(self, directory_name: str):
        '''This function will write a fairly genetic header for the file'''

        with open(self.file_name, "a+") as readme_file:

            readme_file.write(
                f"# README for the {directory_name} directory:\n\n")

    def add_line(self, info_str: str):
        '''This function will the info_str to a new line in the document'''

        with open(
                self.file_name,
                "a+",
        ) as readme_file:
            readme_file.write(info_str + "\n")

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
        logger.info(f"{key} : {value}")
    
