import logging
import getpass
import os
from os import path


def create_logger(name: str, log_filename: str) -> object:
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

    # check if file exists and if it does then deleting the file from the previous run
    if path.exists(log_filename):
        os.remove(log_filename)

    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
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


def record_user_arguments(logger: object, sys_args: object):
    """Function to record the initial arguments provided by the user
    Parameters
    __________
    logger : object
        this is the logger object that information is being written to

    sys_args : object
        this is the object that contains all of the arguments passed by the user

    """

    logger.info(f"Provided Binary file: {sys_args.binary_file}")
    logger.info(f"Recoded options used for PLINK: {sys_args.recode_options}")
    logger.info(
        f"ibd programs used for output within the analysis: {sys_args.ibd_programs}"
    )
    logger.info(
        f"File containing information about the population demographics found at: {sys_args.pop_info}"
    )
    logger.info(
        f"Population code being used within the analysis: {sys_args.pop_code}")
    logger.info(f"output written to: {sys_args.output}")
