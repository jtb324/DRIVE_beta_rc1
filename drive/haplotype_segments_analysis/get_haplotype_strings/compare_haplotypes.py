import pandas as pd
import os

# imports from other files
import population_filter_scripts


def get_haplotype_str(frequency_dict: dict, total_allele_count: int) -> str:
    """function to get the most probable haplotype at each locations
    Parameters
    __________
    frequency_dict : dict
        dictionary containing the number of each allele for each 
        allele in the haplotype from the ped file

    total_allele_count : int
        value of the total alleles present in the dataset for that variant. 
        This is used to normalize the frequencies

    Returns
    _______
    str 
        returns a string containing the most probably haplotype
    """
    probable_haplotype_list: list = []

    haplotype_positions_list: list = list(frequency_dict.keys())

    haplotype_positions_list.sort()

    i, j = 0, 1

    while j <= haplotype_positions_list[-1]:

        #  Creating three list to compare the frequencies within each
        # dictionary and the frequencies across the dictionary

        allele_freq_dict: dict = dict()

        # gets the frequencies for alleles that are in both frequency_dict[i]
        # and frequency_dict[j] as well as the frequencies not in
        # frequency_dict[j]
        for key in frequency_dict[i]:

            if key in frequency_dict[j]:

                allele_freq_dict.setdefault(
                    key, (frequency_dict[i][key] + frequency_dict[j][key]) /
                    (2 * total_allele_count))

            else:

                allele_freq_dict.setdefault(
                    key, frequency_dict[i][key] / (2 * total_allele_count))
        # getting the frequencies that are not in frequency_dict[j]
        for key in frequency_dict[j]:

            if key not in frequency_dict[i]:

                allele_freq_dict.setdefault(
                    key, frequency_dict[j][key] / (2 * total_allele_count))

        probable_haplotype_list.append(
            str(max(allele_freq_dict, key=allele_freq_dict.get)))

        i += 2
        j += 2

    haplotype_str: str = "".join(probable_haplotype_list)

    return haplotype_str


# @population_filter_scripts.pop_filter_decorator
def compare_haplotypes(*args,
                       haplotype: str = None,
                       allele_freq_dict: dict = None) -> str:
    """main function for comparing the frequencies of each allele
    Parameters
    ___________
    file_path : str
        file path to the ped file for either hapibd or ilash that contains 
        the shared segment for all individuals in the mega array data.

    pop_info : str
        file path to a file that list the population code ("EUR", "EAS", etc.)
        for each iid

    pop_code : str
        string that contains the population code to filter the data down from.
        This code should be one of the 1000 genomes code such as "EUR"

    Returns
    _______
    str 
        returns a string of the correct haplotype at each site
    """

    ped_file: str = args[0]

    pop_file_path: str = args[1]

    # load the pop_file into a pandas dataframe
    pop_df: pd.DataFrame = pd.read_csv(pop_file_path, sep="\t")

    pop_code: str = args[2]

    pop_subset_grids: pd.DataFrame = pop_df[pop_df.Pop == pop_code]["grid"]

    with open(ped_file, "r+") as ped_recode_file:

        allele_freq_dict: dict = {}

        # Keeping track of the total number of rows because this is the total
        # number of alleles in the file
        total_alleles: int = 0
        for line in ped_recode_file:

            split_line: list = line.split(" ", 6)

            fid: str = split_line[0]

            haplotype_list: list = split_line[6].split(" ")
            haplotype_len: int = len(haplotype_list)

            if not allele_freq_dict:

                allele_freq_dict = {i: {} for i in range(0, haplotype_len)}

            # check to see if the file is in the specified ancestry
            if not pop_subset_grids.isin([fid]).empty:
                # start a counter
                counter = 0

                for allele in haplotype_list:
                    # inserting the allele into
                    if allele.strip("\n") in allele_freq_dict[counter]:

                        allele_freq_dict[counter][allele.strip("\n")] += 1

                    else:

                        allele_freq_dict[counter].setdefault(
                            allele.strip("\n"), 1)
                    counter += 1

            total_alleles += 1

    # TODO: need function that will compare frequencies of the haplotypes
    haplotype_str: str = get_haplotype_str(allele_freq_dict, total_alleles)

    return haplotype_str
