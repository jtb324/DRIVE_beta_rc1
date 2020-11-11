from os import remove
from os import path
import pandas as pd
import argparse
import os
from pathlib import Path
import shutil
import matplotlib.pyplot as plt
from dataclasses import dataclass

from pandas.core.dtypes.missing import notna


def sort_df(dataframe: pd.DataFrame, column_name: str) -> pd.DataFrame:
    '''This function sort the dataframe based on the provided column name'''

    return dataframe.sort_values(by=[column_name])


def get_unique_variants(dataframe: pd.DataFrame) -> list:
    '''This function will get the unique variants'''

    uniq_var_set: set = set(dataframe.variant_id.tolist())

    return list(uniq_var_set)


def subset_df(dataframe: pd.DataFrame, variant: str) -> pd.DataFrame:
    '''This function will subset the provided dataframe for a specific variant'''

    df_subset: pd.DataFrame = dataframe[dataframe.variant_id == variant]

    return df_subset


def get_len_values(subset_df: pd.DataFrame, column_name: str) -> list:
    '''This function gets a list of all the lengths for each variant'''

    return subset_df[column_name].tolist()


def create_output_dir(output):
    '''This function will create the output directory for the histograms'''

    try:
        os.mkdir(output)

    except FileExistsError:
        pass


def remove_output_dir(output):
    '''This function will remove the output directory from a previous run'''

    output_path = Path(output)

    if output_path.exists() and output_path.is_dir():
        shutil.rmtree(output_path)


def remove_file(file_path: str):
    '''This function will be used to remove files from previous runs'''

    if path.exists(file_path):

        os.remove(file_path)


@dataclass
class length_class:
    program_len_list: list
    list_len: int
    has_NA: bool

    def __init__(self, len_list: list) -> None:

        self.program_len_list = len_list
        self.list_len = len(self.program_len_list)
        self.has_NA = self.is_na()

        if self.has_NA:

            self.program_len_list = [
                element for element in self.program_len_list if element != "N/A"]

            self.list_len = len(self.program_len_list)

    def is_na(self) -> bool:

        is_na_bool: bool = "N/A" in self.program_len_list

        return is_na_bool


def create_histogram(ilash_object, hapibd_object, output_path: str, file_id: str):
    '''This function will create a histogram'''
    print(file_id)
    # creating the file name
    output: str = "".join(
        [output_path, file_id, "_", "segment_distribution", ".png"])

    # creating a subplot

    figure, axes = plt.subplots(nrows=2, ncols=1, sharey="row")

    # setting padding
    figure.tight_layout(pad=3, h_pad=4)

    # creating two values for the axes of each figure
    axes1, axes2 = axes.flatten()

    # making the figure titles

    axes1.set_title("ilash shared segments")
    axes2.set_title("hapibd shared segments")

    # setting the axes labels
    axes1.set_xlabel("Segment Length (cM)")

    axes2.set_xlabel("Segment Length (cM)")

    axes1.set_ylabel("frequency")

    axes2.set_ylabel("frequency")

    # plotting the ilash_object and the hapibd

    axes1.hist(x=ilash_object.program_len_list, bins="auto", alpha=0.7)

    axes2.hist(x=hapibd_object.program_len_list, bins="auto", alpha=0.7)

    figure.savefig(output)

    plt.close()


def remove_na(dataframe: pd.DataFrame) -> pd.DataFrame:
    '''This function will drop all rows that have NaN in the hapibd and ilash len columns'''

    return dataframe[(dataframe["hapibd_len"].notna()) & (dataframe["ilash_len"].notna())]


def write_failed_var(output: str, var: str, chr_num: str):
    '''This function will create a file that list which variants where not able to make a 
    histogram for due to there being no reported ilash or hapibd length'''

    output_full: str = "".join(
        [output, "variants_failed_segment_len_analysis.txt"])

    with open(output_full, "a+") as error_file:

        if os.path.getsize(output) == 0:

            error_file.write(f"variant_id\tchr\n")

        error_file.write(f"{var}\t{chr_num}")


def write_var_info(output: str, var: str, chr_num: str, dataframe_len: str, filtered_df_len: str):
    '''This function will write information about each variant so that I know how many pairs are actually in the run'''
    output_full: str = "".join(
        [output, "histograms_info.txt"])

    with open(output_full, "a+") as info_file:

        info_file.write(
            f"Creating a histogram for variant {var} on chromosome {chr_num}\n")
        info_file.write(f"There were {dataframe_len} pairs identified\n")
        info_file.write(
            f"After removing rows where the haplotype analysis failed, there were {filtered_df_len} pairs left\n")


def run(args):
    "function to run"
    plt.style.use('ggplot')

    remove_output_dir(args.output)

    create_output_dir(args.output)

    # removing the previous files that are formed in anything fails
    remove_file(
        "".join([args.output, "variants_failed_segment_len_analysis.txt"]))

    haplotype_df: pd.DataFrame = pd.read_csv(args.input, sep="\t")

    # sorting the df
    sorted_df: pd.DataFrame = sort_df(haplotype_df, "variant_id")

    uniq_var: list = get_unique_variants(sorted_df)

    # This for loop will subset the dataframe for the specific variant
    for var in uniq_var:

        df_subset: pd.DataFrame = subset_df(sorted_df, var)

        chr_num: str = list(set(df_subset.chr))[0]

        # need to drop rows where both sets are N/A
        filtered_df: pd.DataFrame = remove_na(df_subset)

        # writing some information about how many pairs were identified for the variant
        write_var_info(args.output, var, chr_num, str(
            len(df_subset)), str(len(filtered_df)))

        if filtered_df.empty:

            write_failed_var(args.output, var, chr_num)

            continue

        # This gets a list of lengths for the shared segments from the dataframe
        hapibd_len_object: dataclass = length_class(
            get_len_values(df_subset, "hapibd_len"))
        ilash_len_object: dataclass = length_class(
            get_len_values(df_subset, "ilash_len"))

        # making the hapibd_len_list and ilash_len_list into plots
        # creating a subset directory for the histograms based of of individual variants
        create_output_dir("".join([args.output, "variants/"]))

        create_histogram(ilash_len_object, hapibd_len_object,
                         "".join([args.output, "variants/"]), var)

    # getting the segment lengths for the entire dataframe
    total_hapibd_len = length_class(get_len_values(sorted_df, "hapibd_len"))

    total_ilash_len = length_class(get_len_values(sorted_df, "ilash_len"))

    create_histogram(total_ilash_len, total_hapibd_len,
                     args.output, "all_variants")


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("--input", help="This argument takes the path of the input haplotype_info.txt file",
                        dest="input", type=str, required=True)

    parser.add_argument("--output", help="This argument takes output path",
                        dest="output", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
