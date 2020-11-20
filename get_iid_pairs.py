import pandas as pd
import numpy as np

# getting the individuals from the haplotype_info.txt file


def get_pair_ids(dataframe: pd.DataFrame) -> tuple:
    '''These function will get a list of each pair and return that list'''

    pair_1_list: list = dataframe.pair_1.tolist()

    pair_2_list: list = dataframe.pair_2.tolist()

    chr_num_list: list = dataframe.chr.tolist()

    hapibd_start_list: list = dataframe.hapibd_start.tolist()

    hapibd_end_list: list = dataframe.hapibd_end.tolist()

    ilash_start_list: list = dataframe.ilash_start.tolist()

    ilash_end_list: list = dataframe.ilash_end.tolist()

    hapibd_pair_list: list = [pair for pair in zip(
        pair_1_list, pair_2_list, chr_num_list, hapibd_start_list, hapibd_end_list)]

    ilash_pair_list: list = [pair for pair in zip(
        pair_1_list, pair_2_list, chr_num_list, ilash_start_list, ilash_end_list)]

    return hapibd_pair_list, ilash_pair_list
