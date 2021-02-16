import pandas as pd
import os

# imports from other files
import population_filter_scripts


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

    print("comparing haplotypes")

    ped_file: str = args[0]

    pop_file_path: str = args[1]

    # load the pop_file into a pandas dataframe
    pop_df: pd.DataFrame = pd.read_csv(pop_file_path, sep="\t")

    pop_code: str = args[2]

    pop_subset_grids: pd.DataFrame = pop_df[pop_df.Pop == pop_code]["grid"]

    with open(ped_file, "r+") as ped_recode_file:

        allele_freq_dict: dict = {}

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

                print(allele_freq_dict)
