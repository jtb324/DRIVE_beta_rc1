import sys
import logging
import os
from os import path

# This script will ask for the user to input certain initial parameters


def get_dict_of_variables(class_object: object) -> dict:
    """function to get a list of all the attributes of a class objects
    Parameters
    __________
    class_object: object
        object that is the Input_Gather class

    Returns
    _______
    dict
        returns a dictionary containing all the attributes of the class 
        object
    """
    return class_object.__dict__


class Input_Gather:
    def __init__(self):
        self.logger = self.get_logger()

        self.ANALYSIS_TYPE = None
        self.ask_for_analysis_type()

        self.MIN_CM = None
        self.ask_for_min_cm()

        self.THREADS = None
        self.ask_for_thread_count()

        self.MAF_THRESHOLD = None
        self.ask_for_maf_filter()

    @staticmethod
    def get_logger() -> object:
        """function to get the logger that has the name __main__
        Returns
        _______
        object
            returns the logger object that has the name __main__
        """
        logger = logging.getLogger("__main__")
        return logger

    # Function for analysis type
    def ask_for_analysis_type(self) -> str:
        """
        This function will get the analysis type for the run.
        This can be either "maf", "gene", or ""
        """
        # getting the logger

        ANALYSIS_TYPE: str = input("Please input an analysis type: ").strip(
            " ")

        if ANALYSIS_TYPE != "":
            self.logger.info(f"Using the analysis type: {ANALYSIS_TYPE}")
        else:
            self.logger.info("no analysis type passed to the program")

        if ANALYSIS_TYPE not in ["gene", "maf", ""]:

            print(
                "Please choose one of the allowed analysis types: 'gene'/'maf'/'' "
            )

            self.logger.error(
                "invalid analysis type passed to the program. Terminating program"
            )
            sys.exit(1)

        self.ANALYSIS_TYPE = ANALYSIS_TYPE

    def ask_for_min_cm(self) -> int:
        """
        This function will ask the user for the minimum centimorgan
        threshold to be used in analysis
        """
        MIN_CM: str = input(
            "Please input a value for the minimum centimorgan threshold: (Default is 3 cM) "
        )

        if not MIN_CM:
            MIN_CM = 3
        else:
            try:
                MIN_CM = int(MIN_CM)
            except ValueError:
                print(
                    "invalid parameter passed for the minimum centimorgan threshold."
                )

                self.logger.error(
                    "Invalid parameter passed for the minimum centimorgan threshold. Please input an integer value"
                )

        self.logger.info(f"Minimum Centimorgan Threshold: {MIN_CM}")

        self.MIN_CM = MIN_CM

    def ask_for_thread_count(self) -> int:
        """
        This function will ask the user for the number of threads to be used 
        during the analysis. It returns an integer. The function defaults to 
        3 if the user provides no input
        """
        THREADS: str = input(
            "Please enter the number of threads you wish to use during this process. The default value is 3. (Bear in mind that this number will be used for all parallelized steps): "
        )
        # creating a default value
        if not THREADS:

            THREADS = 3

        else:
            # checking to see if the provided value can be converted to an integer value
            try:
                THREADS = int(THREADS)
            except ValueError:
                print(
                    "invalid parameter passed for the number of threads to be used during the analysis"
                )
                self.logger.error(
                    "invalid parameter passed for the number of threads to be used during the analysis. Please enter an number"
                )

        self.logger.info(
            f"Number of Threads used during the computation: {THREADS}")

        self.THREADS = THREADS

    def ask_for_maf_filter(self) -> str:
        """
        This function will ask for a minor allele frequency threshold from the 
        user. If the user provides no value then it will default to 0.05. The
        function returns a string
        """
        MAF_FILTER: str = input(
            "Please input a minor allele frequency threshold to be used. (The default is 0.05): "
        )

        if not MAF_FILTER:

            MAF_FILTER: str = '0.05'
        if float(MAF_FILTER) > 0.5:
            print(
                "minor allele frequency filter must be lower than 0.5. System terminating...")

            self.logger.error(
                f"A minor allele frequency of {MAF_FILTER} was passed to the program minor allele frequency filter must be lower than 0.5.")

            sys.exit(1)

        self.logger.info(f"MAF_FILTER: {MAF_FILTER}")

        self.MAF_THRESHOLD = MAF_FILTER
