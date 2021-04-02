import sys
import pandas as pd
import os


def load_frequency_file(file_path: str) -> pd.DataFrame:
    '''This function loads the provided file into a pandas dataframe'''
    try:
        return pd.read_csv(file_path, sep="\t")
    except FileNotFoundError:
        raise Exception("file was not found")


def filter_for_higher_maf(frequency_df: pd.DataFrame,
                          threshold: int) -> pd.DataFrame:
    '''This function filters the provided dataframe for values where the
    allele frequency is greater than the threshold'''

    return frequency_df[frequency_df.allele_freq > threshold]


def get_filtered_var_list(df_subset: pd.DataFrame) -> list:
    '''This function will return a list of all the variants that exceeded 
    the threshold'''

    return df_subset.variant_id.values.tolist()


def check_mafs(file_path: str, threshold: int) -> tuple:
    '''This function will check to see if any variants exceed a threshold
    the provided population file'''

    frequencies_df: pd.DataFrame = load_frequency_file(file_path)

    df_subset: pd.DataFrame = filter_for_higher_maf(frequencies_df, threshold)

    variant_list: list = get_filtered_var_list(df_subset)

    escape_char_pressed: bool = False

    program_end_code: int = 1
    # if the variant_list is not empty
    if variant_list:

        print(
            f"The variants, {', '.join(variant_list)}, exceed the provided threshold of {threshold}"
        )
        print(
            "It is recommended to remove these variants. The user can leave these variants in but the programs runtime may be longer than expected."
        )
        while not escape_char_pressed:

            user_input: str = input("Would you like to continue? y/n: ")

            if user_input == "y":

                program_end_code = 1

                escape_char_pressed = True
            elif user_input == "n":

                program_end_code = 0

                escape_char_pressed = True
            else:
                print("invalid option selected")

        return (variant_list, program_end_code)
    else:
        return (variant_list, program_end_code)
