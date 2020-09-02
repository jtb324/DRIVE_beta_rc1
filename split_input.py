# This file is used if an input file containing multiple chromosomes is used.
# This file will just split the file into different list of variants based off of
# the chromosomes

# Need to have python package "xlrd" installed. Just run conda/pip install xlrd

# Import packages
import pandas as pd
import argparse
import sys

from check_directory import check_dir
from file_exist_checker import Check_File_Exist
from plink_python_run import PLINK_Runner


class Input_Chr_Splitter(Check_File_Exist):
    '''This class splits the input file which contains variants for multiple chromosomes into
    files that contain variants for just one chromosome'''

    def __init__(self, variant_csv_path: str, output_path: str):

        # Create property called file that is just the path to the input file
        self.file = variant_csv_path

        # check if the file is an csv file
        if self.file[-4:] == ".csv":

            self.var_df = self.check_file_exist(separator=",")

        # check if it is an excel file
        elif self.file[-5:] == ".xlsx":

            self.var_df = pd.read_excel(self.file)

        # Fail if it is any other file format
        else:
            print(f"The file {self.file} is not a supported file type. \
                Supported file types are .xlsx and .csv")

        self.output_path = check_dir(output_path, "variants_of_interest/")

        self.split_input_file()

    def split_input_file(self):

        # Get a list of all chromosome numbers in file

        chr_list = self.var_df.Chr.unique().tolist()

        # Iterate through chr_list to isolate dataframe for a specific chromosome
        for chromo in chr_list:

            variant_df_subset = self.var_df[self.var_df["Chr"] == chromo]

            variant_df_subset = variant_df_subset.reset_index(drop=True)

            # writing this  subset to csv file

            variant_df_subset.to_csv(
                "".join([self.output_path, "variants_of_interest", ".chr", str(chromo), "_list", ".csv"]))

            self.write_variants_to_file(variant_df_subset, str(chromo))

    def write_variants_to_file(self, variant_df_subset: pd.DataFrame, chromosome: str):
        # Now create a list of just the variants for that chromosome
        # Need to isolate the SNP column for the subset df
        variant_list = variant_df_subset.SNP.values.tolist()

        # write the variant_list to a file
        # Open a file at the out
        MyFile = open(
            "".join([self.output_path, "variants_of_interest", ".chr", chromosome, "_list", ".txt"]), 'w')

        for variant_id in variant_list:
            # Write each SNP id to each line in the txt file
            MyFile.write(variant_id)
            MyFile.write('\n')
        MyFile.close()


def run(args):
    print("splitting the single file of multiple chromosomes into multiple files of a single chromosome")

    Input_Chr_Splitter(args.input, args.output)

    print("running PLINK...")

    for recode_option in args.recode:
        recode_flag = "".join(["--", recode_option])
        plink_runner = PLINK_Runner(args.bfile, args.var_dir, recode_flag)


def main():
    parser = argparse.ArgumentParser(
        description="This CLI splits the input file into separate files for each chromosome")

    parser.add_argument("--input", help="This argument takes the initial input file which contains variants for multiple chromosomes and splits it into multiple files, one for each chromosome",
                        dest="input", type=str, required=True)

    parser.add_argument("--output", help="This argument creates the main output directory for the chromosome files",
                        dest="output", type=str, required=True)

    parser.add_argument("--bfile", help="This argument leads PLINK to the correct binary file path",
                        dest="bfile", type=str, required=True)

    parser.add_argument("--var_list_dir", help="This argument provides the directory of all the files that lis the variants",
                        dest="var_dir", type=str, required=True)

    parser.add_argument("--recode", help="This argument tells what kind of recode you want PLINK to run",
                        dest="recode", type=str, nargs="+", required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
