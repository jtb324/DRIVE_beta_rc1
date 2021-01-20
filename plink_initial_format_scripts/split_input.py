# This file is used if an input file containing multiple chromosomes is used.
# This file will just split the file into different list of variants based off of
# the chromosomes

# Need to have python package "xlrd" installed. Just run conda/pip install xlrd

# Import packages
import pandas as pd
import argparse
import sys

import file_creator_scripts
import plink_initial_format_scripts


class Input_Chr_Splitter:
    """This class splits the input file which contains variants for multiple chromosomes into
    files that contain variants for just one chromosome"""
    def __init__(self, variant_csv_path: str, output_path: str):

        # Create property called file that is just the path to the input file

        self.file = variant_csv_path
        # check if the file is an csv file
        if self.file[-4:] == ".csv":

            self.var_df = pd.read_csv(self.file)

        # check if it is an excel file
        elif self.file[-5:] == ".xlsx":

            self.var_df = pd.read_excel(self.file)

        # Fail if it is any other file format
        else:
            print(f"The file {self.file} is not a supported file type. \
                Supported file types are .xlsx and .csv")

        self.output_path = file_creator_scripts.check_dir(
            output_path, "plink_output_files/")

        self.split_input_file()

    def split_input_file(self):

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
                self.output_path,
                "variants_of_interest",
                ".chr",
                chromo,
                "_list",
                ".csv",
            ]))

            self.write_variants_to_file(variant_df_subset, str(chromo))

    def write_variants_to_file(self, variant_df_subset: pd.DataFrame,
                               chromosome: str):
        # Now create a list of just the variants for that chromosome
        # Need to isolate the SNP column for the subset df
        variant_list = variant_df_subset.SNP.values.tolist()

        # write the variant_list to a file
        # Open a file at the out
        MyFile = open(
            "".join([
                self.output_path,
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


def split_input_and_run_plink(
    input_file: str,
    output: str,
    recode_options: list,
    binary_file: str,
    plink_files_dir: str,
) -> str:
    print(
        "splitting the single file of multiple chromosomes into multiple files of a single chromosome"
    )

    Input_Chr_Splitter(input_file, output)

    print("running PLINK...")

    plink_runner = plink_initial_format_scripts.PLINK_Runner(
        recode_options, output, binary_file, var_list_dir=plink_files_dir)

    variant_file_list: list = plink_runner.generate_file_list()

    return plink_runner.run_PLINK_snps(variant_file_list)
