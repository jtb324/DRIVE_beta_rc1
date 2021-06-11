import sys
import pandas as pd
import os
from typing import Dict, Tuple, List
from functools import partial


def check_frequencies(threshold: int, variant_freq: Tuple) -> bool:
    """function to check if a variant allele frequency exceeds the provided threshold
    Parameters
    __________
    threshold : int
        integer value for the threshold. This is provided intially by 
        the user
    
    variant_freq : Tuple
        tuple where the first value is the variant id and the second 
        value is the frequency

    Returns
    _______
    bool
        returns true or false based on whether the variant is <= to the 
        threshold
    """
    freq: float = variant_freq[1]

    return freq > threshold

def check_mafs(maf_dict: Dict, threshold: int) -> tuple:
    '''This function will check to see if any variants exceed a threshold
    the provided population file'''

    # creating a list to keep all of the variant names above the 
    # threshold
    variant_list: List[str] = []

    # iterating through all the variants to determine which of the 
    # allele frequencies are higher than the threshold. These variants 
    # are added to teh variant_list function
    for _, variant_maf_dict in maf_dict.items():

        var_above_threshold: List[Tuple] = list(filter(partial(check_frequencies, threshold), variant_maf_dict))

        var_names: List[str] = [var_tuple[0] for var_tuple in var_above_threshold] 

        variant_list = variant_list + var_names


    escape_char_pressed: bool = False

    program_end_code: int = 1

    # if the variant_list is not empty then the program suggest that 
    # the user pause the program and remove these. If the user wants to
    # continue then they can user can enter y but the runtime my 
    # increase
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
