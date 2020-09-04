# This file is used to run plink from python

import subprocess
import pandas as pd
import os
import sys
import glob


class PLINK_Runner:
    def __init__(self, binary_file: str, var_list_directory: str, recode_flag: str):
        self.binary_file = binary_file
        self.current_dir = os.getcwd()
        self.var_list_dir = var_list_directory
        self.recode = recode_flag

        self.change_directory()

        var_file_list = self.generate_file_list()

        self.run_PLINK(var_file_list, self.recode)

    def change_directory(self):
        '''This function switches to the directory of the variants specified'''

        os.chdir(self.var_list_dir)

    def generate_file_list(self) -> list:
        '''This function will return a list of all the variant files that can be fed to PLINK'''

        file_list = []

        for file in glob.glob("*.txt"):

            if len(file) == 0:
                print(
                    "There were no txt files found which contained a list of variant ids to be fed to PLINK")

                sys.exit(1)

            full_file_path = "".join([self.var_list_dir, file])

            print(full_file_path)

            file_list.append(full_file_path)

        return file_list

    def run_PLINK(self, file_list: list, recode_option: str):
        '''This function will use the subprocess module to run PLINK'''

        os.chdir(self.current_dir)

        for var_file in file_list:

            output_file_name = var_file[:-4]

            subprocess.run(["plink",
                            "--bfile",
                            self.binary_file,
                            "--extract",
                            var_file,
                            "--out",
                            output_file_name,
                            recode_option,
                            ])
