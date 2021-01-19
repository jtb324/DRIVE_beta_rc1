# This file is used to run plink from python

import subprocess
import os
import sys
import glob


class PLINK_Runner:
    def __init__(
        self,
        recode_flag: list,
        output: str,
        binary_file: str,
        **name,e
    ):
        self.binary_file = binary_file
        self.output = output
        self.current_dir = os.getcwd()
        self.var_list_dir = name["var_list_dir"]
        self.recode = recode_flag
        if "maf_filter" in name:

            self.maf = name["maf_filter"]
        if "start" in name:

            self.start = name["start"]
        if "end" in name:

            self.end = name["end"]

    def generate_file_list(self) -> list:
        """This function will return a list of all the variant files that can be fed to PLINK"""

        os.chdir(self.var_list_dir)

        file_list = []

        for file in glob.glob("*.txt"):

            if len(file) == 0:
                print(
                    "There were no txt files found which contained a list of variant ids to be fed to PLINK"
                )

                sys.exit(1)

            full_file_path = "".join([self.var_list_dir, file])

            file_list.append(full_file_path)

        os.chdir(self.current_dir)

        return file_list

    def run_PLINK_snps(self, file_list: list):
        """This function will use the subprocess module to run PLINK and extract snps from a specified list"""

        for var_file in file_list:

            for option in self.recode:
                print(option)
                output_file_name = var_file[:-4]
                subprocess.run(
                    [
                        "plink",
                        "--bfile",
                        self.binary_file,
                        "--extract",
                        var_file,
                        "--out",
                        output_file_name,
                        "".join(["--", option]),
                    ],
                    check=False,
                )

        return "".join([self.output, "plink_output_files/"])

    def run_PLINK_maf_filter(self) -> str:
        """This function will use the subprocess module to run PLINK and extract snps from a specified list"""

        full_output_path: str = "".join(
            [self.output, "plink_output_files/", self.start, "_", self.end])

        for options in self.recode:
            if self.start and self.end:
                subprocess.run(
                    [
                        "plink",
                        "--bfile",
                        self.binary_file,
                        "--max-maf",
                        self.maf,
                        "--out",
                        full_output_path,
                        "--from",
                        self.start,
                        "to_rs",
                        self.end,
                        "".join(["--", options]),
                    ],
                    check=False,
                )
            else:
                subprocess.run(
                    [
                        "plink",
                        "--bfile",
                        self.binary_file,
                        "--max_maf",
                        self.maf,
                        "--out",
                        self.output,
                        "".join(["--", options]),
                    ],
                    check=False,
                )
        return "".join([self.output, "plink_output_files/"])
