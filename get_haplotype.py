import pandas as pd
import glob
import os
import argparse


def get_files(file_dir: str, file_extension) -> list:
    '''This function get the files for the directory'''

    cur_dir: str = os.getcwd()

    os.chdir(file_dir)

    file_list = []

    for file in glob.glob(file_extension):

        full_file_path = "".join([file_dir, file])

        file_list.append(full_file_path)

    os.chdir(cur_dir)

    return file_list


class carrier_files:
    '''This class loads the pair list and gets the unique variants from it'''

    def __init__(self, pair_file: str):
        self.pairs_file = pair_file

    def load_confirmed_carrier_files(self) -> pd.DataFrame:
        '''This function will load the provided file into a dataframe'''

        carriers_df: pd.DataFrame = pd.read_csv(self.pairs_file, sep="\t")

        return carriers_df

    @staticmethod
    def get_unique_variants(pairs_df: pd.DataFrame) -> list:
        '''This function will return a list of all the unique variants'''

        # the following line gets a list of all the unique variants within the dataframe
        variant_list: list = pairs_df.variant_id.values.unique().tolist()

        return variant_list

    @staticmethod
    def get_unique_network_ids(carriers_df_subset: pd.DataFrame) -> list:
        '''This function will take a subset of the carriers_df only for a specific variant and
        return a list of all the unique network ids'''

        # the following line gets a list of all the unique variants within the dataframe
        id_list: list = carriers_df_subset["Network ID"].values.unique(
        ).tolist()

        return id_list


def get_index_positions(file_path: str) -> list:

    # This dictionary contains the index of values such as IID1, IID2, Start point
    # end point, and the shared segment length. True will contain indixes for an iLash
    # file and false will contain the indices for a hapibd file
    indx_dict: dict = {
        True: [0, 2, 5, 6, 9],
        False: [0, 2, 5, 6, 7]

    }

    indx_list: list = indx_dict["match" in file_path]

    return indx_list


def get_pair_list(df_subset: pd.DataFrame, col_num: int) -> list:
    '''This function will take the df subset and and specific column such as 1 or 2 and return 
    a list of all the pairs for that column'''

    col_handler: dict = {
        1: "pair_1",
        2: "pair_2"
    }

    # This next line gets the name of the column
    col_name: str = col_handler[col_num]

    # This next line gets a list of either all the pair 1 or pair 2 values
    pair_list: list = df_subset[col_name].values.tolist()

    return pair_list


def get_chr_num(dataframe: pd.DataFrame) -> str:
    '''This function gets the chromosome number for the specific variant and network id'''

    # This gets the unique values of the chromosome numbers
    chr_num: str = str(dataframe.chr_num.values.unique().tolist())

    return chr_num


def get_ibd_file(file_list: list, chr_str: str) -> str:
    '''This function gets the ibd files using list comprehension'''

    ibd_file: str = [file for file in file_list if chr_str in file][0]

    return ibd_file


def get_haplotype_info(ibd_file: str, pair_1: str, pair_2: str) -> dict:
    '''This function will get the haplotype information for the pair in the ibd file'''

    indx_list: list = get_index_positions(ibd_file)

    # Getting the indices for the positions of interest
    pair_1_indx: int = indx_list[0]
    pair_2_indx: int = indx_list[1]
    start_indx: int = indx_list[2]
    end_indx: int = indx_list[3]
    segment_length_indx: int = indx_list[0]


def run(args):
    "function to run"
    # getting the suffix for ilash files from one of the arguments passed
    ilash_file_suffix: str = [
        ending for ending in args.suffix if "match" in ending][0]

    # getting the suffix for hapibd files from one the arguments passed
    hapibd_file_suffix: str = [
        ending for ending in args.suffix if "ibd" in ending][0]

    # getting the list of ilash files
    ilash_file_list: list = get_files(args.ilash, ilash_file_suffix)

    # getting the list of hapibd files
    hapibd_file_list: list = get_files(args.hapibd, hapibd_file_suffix)

    carrier_file_handler = carrier_files(args.carriers_file)

    carriers_df: pd.DataFrame = carrier_file_handler.load_confirmed_carrier_files()

    # Getting a list of all the unique variants
    unique_var_list: list = carrier_file_handler.get_unique_variants(
        carriers_df)

    # iterating through each variant
    for variant in unique_var_list:

        # Restricting the carrier dataframe to only those variants
        carrier_df_subset: pd.DataFrame = carriers_df[carriers_df.variant_id == variant]

        # Getting a list of all the network ids
        unique_network_ids: list = carrier_files.get_unique_network_ids(
            carrier_df_subset)

        for network_id in unique_network_ids:

            # subsetting the dataframe further just for individuals with a certain network id
            in_network_df: pd.DataFrame = carrier_df_subset[carrier_df_subset["Network ID"] == network_id]

            # getting all of the first pairs in a list
            pair_1_list: list = get_pair_list(in_network_df, 1)

            # getting all of the second pairs in a list
            pair_2_list: list = get_pair_list(in_network_df, 2)

            # getting the chromosome number
            chr_num: str = get_chr_num(in_network_df)

            # adding chr and . to the chr_num so that it matches the files
            chr_str: str = "".join(["_chr", chr_num, "."])

            # getting the specific ilash file for that chromosome
            ilash_ibd_file: str = get_ibd_file(ilash_file_list, chr_str)

            # getting the specific hapibd file for that chromosome
            hapibd_ibd_file: str = get_ibd_file(hapibd_file_list, chr_str)

            # iterating through the pairs list
            for pair_1, pair_2 in zip(pair_1_list, pair_2_list):


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("--ilash", help="This argument list the path to the directory of ilash files",
                        dest="ilash", type=str, required=True)

    parser.add_argument("--hapibd", help="This argument list the path to the directory of hapibd files",
                        dest="hapibd", type=str, required=True)

    parser.add_argument("--suffix", help="This argument takes a list of the possible suffixes of each file output from the ibd files",
                        dest="suffix", type=str, nargs="+", required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
