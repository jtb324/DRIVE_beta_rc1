# This file will contain the parent class for generating readme and log files
import os
from datetime import datetime


class Documentation:
    '''This class will be the parent class for the documentation throughout the program'''

    def __init__(self, file_name: str, output_path: str):

        self.file_name = "".join([output_path, file_name])

    def create_file(self):
        with open(self.file_name, "w+") as file:
