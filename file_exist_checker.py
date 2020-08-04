################################################
# Import modules
import sys
from os import path
import pandas as pd

################################################


class Check_File_Exist:

    def __init__(self, file_to_check, logger):
        '''This function test to see if the passed file exist'''

        self.file = file_to_check
        self.log_file = logger

    def check_file_exist(self, column_names=None, separator=" ", header_value='infer', skip_rows=0):
        '''This function checks if the file exist and then will either return the loaded file if it does exist or it will return an error message saying the file was note found.'''

        try:

            file = pd.read_csv(self.file, sep=separator,
                               names=column_names, header=header_value, skiprows=skip_rows)

        except FileNotFoundError:

            print("The raw recoded file at {} was not found.".format(
                self.file))

            self.log_file.error("The file at {} was not found.".format(
                self.file))

            sys.exit(1)

        return file
