# This script is designed to make readmes for the directory that it is in

import os
from os import path


class Readme_Generator():
    '''This class will generate READMEs for subdirectories'''
    def __init__(self, output_dir: str, text: str, header: str, **name):

        self.output_path: str = "".join([output_dir, "_README.md"])
        self.header: str = header
        self.text: str = text

        if "binary_file" in name:
            self.binary = name["binary_file"]

    def check_if_file_exist(self):
        '''Checks if the README already exists'''
        if path.exists(self.output_path):
            os.remove(self.output_path)

    def create_readme(self):
        '''Creates teh readme file'''
        with open(self.output_path, "w+") as readme:

            readme.write(self.text)
