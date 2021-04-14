import os
from os import path
import logging
import sys
import glob
import subprocess
import pandas as pd

from .check_missing_variants import check_for_missing_var
import utility_scripts


class Analysis_Checker:
    def __init__(self, analysis_type: str, recode_options: list,
                 binary_file: str, output_path: str, plink_dir: str, **name):
        self.analysis_type = analysis_type
        self.logger = self.initialize_logger()
        self.recode_flags: list = recode_options
        self.binary_file: str = binary_file
        self.output: str = output_path
        self.plink_dir = plink_dir

        if "var_file" in name:
            self.var_file = name["var_file"]

        if "maf_filter" in name:

            self.maf = name["maf_filter"]
        self.check_if_path_exists()

    def check_if_path_exists(self):
        """function check if the plink directory already exist. If it doesn't then it makes the directory"""

        if not path.exists(self.plink_dir):
            os.mkdir(self.plink_dir)

    @staticmethod
    def initialize_logger() -> object:
        """function that will get the logger with the name __main__
        Returns
        _______
        object
            returns a log object that comments will be logged to
        """
        return logging.getLogger("__main__")

    @utility_scripts.class_readme_generator
    def check_analysis(self, **name):
        """Function to check the analysis type and then run the corresponding plink steps
        """
        # setting initial parameters for the range
        self.START_RS: str = None
        self.END_RS: str = None

        # overwriting the values of START_RS and END_RS if values are present
        if "range" in name:
            self.START_RS = name["range"][0]
            self.END_RS = name["range"][1]

        analysis_handler: dict = {
            "gene": self.gene_analysis,
            "maf": self.maf_function,
            "": self.empty_analysis
        }

        analysis_handler[self.analysis_type]()

    def check_missing_var_count(self) -> int:
        missing_var_count: int = check_for_missing_var(
            self.plink_dir, self.var_file)

        self.logger.info(
            f"There were {missing_var_count} missing variants between the provided file at {self.var_file} and the plink output in the directory {self.plink_dir}"
        )

    def empty_analysis(self):
        self.logger.warning("no analysis type passed to the program")
        self.logger.warning(
            "analysis will assume that the user has provided the necessaryped, map, and raw files from PLINK"
        )

    def maf_function(self):

        # next four lines log information about the directory that the plink
        # output will be written to as well as where the PLINK stdout will
        # be written to
        self.logger.info(
            f"Writting the output of plink to files at {self.plink_dir}")

        self.logger.info(
            f"Writting PLINK's stdout stream to {''.join([self.plink_dir,'plink_log.log'])}"
        )
        self.logger.info(
            f"Beginning analysis for the range starting with the variant {self.START_RS} and ending at {self.END_RS}"
        )

        # getting a chromosome number
        CHR: str = input(
            "Please input the chromosome that the variant of interest is on. Please use a leading 0 for single digit numbers: "
        )
        self.logger.info(
            f"setting the chromosome of interest to be chromosome {CHR}")

        full_output_path: str = "".join([
            self.output, "plink_output_files/", self.START_RS, "_",
            self.END_RS, ".chr", CHR, "_list"
        ])

        with open("".join([self.plink_dir, "plink_log.log"]),
                  "a+") as plink_log:

            for options in self.recode_flags:
                if self.START_RS and self.END_RS:
                    subprocess.run([
                        "plink",
                        "--bfile",
                        self.binary_file,
                        "--max-maf",
                        self.maf,
                        "--out",
                        full_output_path,
                        "--from",
                        self.START_RS,
                        "--to",
                        self.END_RS,
                        "".join(["--", options]),
                    ],
                                   check=False,
                                   stdout=plink_log,
                                   stderr=plink_log)
                else:
                    subprocess.run([
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
                                   stdout=plink_log,
                                   stderr=plink_log)

        self.logger.info(f"PLINK output files written to: {self.plink_dir}")

    # TODO: Create a function to handle the gene analysis

    def gene_analysis(self):
        """This function will use the subprocess module to run PLINK and extract snps from a specified list"""

        self.logger.info(
            f"generating a list of snps from the file: {self.var_file}")

        # creating a object to read the csv file
        Input_handler: object = Input_Splitter(self.var_file, self.output,
                                               self.plink_dir)
        Input_handler.split_input_file(self.plink_dir)
        variant_file_list: list = Input_handler.generate_file_list(
            self.plink_dir)
        for variant_file in variant_file_list:
            with open("".join([self.plink_dir, "plink_log.log"]),
                      "a+") as plink_log:

                for var_file in variant_file_list:

                    for option in self.recode_flags:
                        output_file_name = var_file[:-4]

                        subprocess.run([
                            "plink",
                            "--bfile",
                            self.binary_file,
                            "--max-maf",
                            self.maf,
                            "--extract",
                            var_file,
                            "--out",
                            output_file_name,
                            "".join(["--", option]),
                        ],
                                       check=False,
                                       stdout=plink_log,
                                       stderr=plink_log)

                plink_log.close()
        self.logger.info(f"PLINK output files written to: {self.plink_dir}")


class Input_Splitter:
    """This class splits the input file which contains variants for multiple
    chromosomes into files that contain variants for just one chromosome"""
    def __init__(self, variant_csv_path: str, output_path: str,
                 variant_list_dir: str):
        # Create property called file that is just the path to the input file
        self.logger: object = self.initialize_logger()
        self.cur_dir: str = os.getcwd()
        self.output: str = output_path
        self.file: str = variant_csv_path
        self.var_list_dir: str = variant_list_dir
        # check if the file is an csv file
        if self.file[-4:] == ".csv":

            self.var_df = pd.read_csv(self.file)

        # check if it is an excel file
        elif self.file[-5:] == ".xlsx":

            self.var_df = pd.read_excel(self.file)

            # Fail if it is any other file format and log this failure to a file
        else:
            self.logger.error(
                f"The file {self.file} is not a supported file type.\
                              Supported file types are xlsx and .csv")

            print(f"The file {self.file} is not a supported file type. \
                Supported file types are .xlsx and .csv")

            sys.exit(1)

    @staticmethod
    def initialize_logger() -> object:
        """function that will get the logger with the name __main__
        Returns
        _______
        object
            returns a log object that comments will be logged to
        """
        return logging.getLogger("__main__")

    def split_input_file(self, plink_dir: str):

        # Get a list of all chromosome numbers in file

        chr_list = self.var_df.Chr.unique().tolist()

        # Iterate through chr_list to isolate dataframe for a specific chromosome
        for chromo in chr_list:

            variant_df_subset = self.var_df[self.var_df["Chr"] == chromo]

            variant_df_subset = variant_df_subset.reset_index(drop=True)

            # Making sure that if the chromosome # is a single digit that a leading
            # zero gets added to it
            if len(str(chromo)) == 1:
                chromo = str(chromo).zfill(2)
            else:
                chromo = str(chromo)

            # writing this  subset to csv file
            variant_df_subset.to_csv("".join([
                plink_dir,
                "variants_of_interest",
                ".chr",
                chromo,
                "_list",
                ".csv",
            ]))

            self.write_variants_to_file(variant_df_subset, str(chromo),
                                        plink_dir)

    def write_variants_to_file(self, variant_df_subset: pd.DataFrame,
                               chromosome: str, plink_dir: str):
        # Now create a list of just the variants for that chromosome
        # Need to isolate the SNP column for the subset df
        variant_list = variant_df_subset.SNP.values.tolist()

        # write the variant_list to a file
        # Open a file at the out
        MyFile = open(
            "".join([
                plink_dir,
                "variants_of_interest",
                ".chr",
                chromosome,
                "_list",
                ".txt",
            ]),
            "w",
        )

        for variant_id in variant_list:
            # Write each SNP id to each line in the txt file
            MyFile.write(variant_id)
            MyFile.write("\n")
        MyFile.close()

    def generate_file_list(self, plink_dir: str) -> list:
        """This function will return a list of all the variant files that 
        can be fed to PLINK"""
        # changing to the directory where the variant text files are
        os.chdir(plink_dir)

        file_list = []

        for file in glob.glob("*list.txt"):
            # if no file is present then it will log this as an error
            if len(file) == 0:
                self.logger.info(
                    "There were no txt files found which contained a list of variant ids to be fed to PLINK"
                )

                sys.exit(1)

            file_list.append("".join([plink_dir, file]))

        os.chdir(self.cur_dir)

        return file_list
