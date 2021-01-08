# This file is used to run plink from python

import subprocess
import os
import sys
import glob


class PLINK_Runner:
    def __init__(
        self,
        binary_file: str,
        recode_flag: str,
        output: str,
        maf_filter: str = None,
        var_list_directory: str = None,
    ):
        self.binary_file = binary_file
        self.output = output
        self.current_dir = os.getcwd()
        self.var_list_dir = var_list_directory
        self.recode = recode_flag
        self.maf = maf_filter

    def generate_file_list(self) -> list:
        '''This function will return a list of all the variant files that can be fed to PLINK'''

        os.chdir(self.var_list_dir)

        file_list = []

        for file in glob.glob("*.txt"):

            if len(file) == 0:
                print(
                    "There were no txt files found which contained a list of variant ids to be fed to PLINK"
                )

                sys.exit(1)

            full_file_path = "".join([self.var_list_dir, file])

            print(full_file_path)

            file_list.append(full_file_path)

        os.chdir(self.current_dir)

        return file_list

    def run_PLINK_snps(self, file_list: list):
        '''This function will use the subprocess module to run PLINK and extract snps from a specified list'''

        for var_file in file_list:

            output_file_name = var_file[:-4]

            subprocess.run([
                "plink",
                "--bfile",
                self.binary_file,
                "--extract",
                var_file,
                "--out",
                output_file_name,
                self.recode,
            ])

    def run_PLINK_maf_filter(self):
        '''This function will use the subprocess module to run PLINK and extract snps from a specified list'''

        subprocess.run([
            "plink",
            "--bfile",
            self.binary_file,
            "--max_maf",
            self.maf,
            "--out",
            self.output,
            self.recode,
        ])
