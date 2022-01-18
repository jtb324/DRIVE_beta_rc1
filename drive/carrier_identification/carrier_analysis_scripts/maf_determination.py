import pandas as pd
import glob
import os
from os import path
import re
from typing import Dict, List
import population_filter_scripts
import utility_scripts
# need to gather all of the single var list

def get_allele_frq(carrier_dict: Dict, raw_dir: str, pop_info_file: str, pop_code: str, output: str) -> None:
    """Function to determine the minor allele frequency within the MEGA array dataset
    Parameters
    __________
    carrier_dict : Dict
        Dictionary of dictionaries that has the carriers for each variant for each chromosome

    raw_dir : str
        string that list the directory that the raw files from PLINK 
        are in

    pop_info_file : str
        file that contains information about the ethinicity of each GRID ID within the dataset
    
    pop_code : str
        string that tells what ethnicity the user wishes to filter for

    output : str
        string that list the root output directory for the program to output to

    Returns
    _______
    Dict[str, Dict[str, float]]
        Dictionary that has the minor allele frequency for each variant for each chromosome
    """
    raw_file_list: List[str] = utility_scripts.get_file_list(raw_dir, "*.raw")

    # checking if the output file does exist and removing it if it 
    # does
    utility_scripts.check_file(os.path.join(output, "carrier_analysis_output/allele_frequencies.txt"))

    # opening the file to write to it
    with open(os.path.join(
        output, "carrier_analysis_output/allele_frequencies.txt"), "a+") as myFile:

        # writing a header line
        myFile.write("chr\tvariant_id\tallele_freq\n")

        # creating a dictionary that will be used as a datastructure to 
        # return the maf for each variant
        maf_dict: Dict[str, Dict[str, float]] = {}

        for chr_num, variant_dict in carrier_dict.items():
            
            # adding the chromosome number into the maf_dict
            maf_dict.setdefault(chr_num, {})

            # adding a _ to the chromosome number to match the raw_file
            chr_num = "".join([chr_num, "_"])

            # selecting the correct raw file based on the chromosome 
            # number
            raw_file: str = [
                raw_file for raw_file in raw_file_list
                if chr_num in raw_file
            ][0]

            # loading in the raw file into a dataframe and filtering it just for the desired population code using the Pop_Filter code
            raw_file_df: pd.DataFrame = population_filter_scripts.run_pop_filter(pop_info_file, raw_file,
                                                    pop_code)
            
            # iterating through all the variants
            for variant, iid_list in variant_dict.items():
            
                # determining the total number of alleles by just 
                # multipling the number of rows by two
                total_allele_count: int = raw_file_df.shape[0] * 2

                # subsetting the raw_file_df for the carriers for a specific variant
                carry_allele_df: pd.Series = raw_file_df[raw_file_df.   IID.isin(iid_list)][variant]

                # all of the iids are carriers. This line determines 
                # how many are coded as 1 and how many are coded as 2
                minor_allele_counts: pd.Series = carry_allele_df.value_counts()

                # setting an initial counter for the minor_allele_count
                minor_allele_count: int = 0
                
                # iterating through the value counts which will be
                # either 1.0 or 2.0. This allows the minor allele count
                # to be determined
                for index, counts in minor_allele_counts.iteritems():
                    
                    # creating a handler that will either return the 
                    # counts or double the allele counts value if the 
                    # the index represents homozygous recessive counts
                    index_handler: Dict = {
                        True : counts,
                        False : 2 * counts
                    }
                    
                    minor_allele_count += index_handler[index == 1.0] 
                    
                minor_allele_frq: float = minor_allele_count / total_allele_count

                myFile.write(
                    f"{chr_num[:len(chr_num)-1]}\t{variant}\t{minor_allele_frq}\n")

                # adding the allele frequency for each variant into the 
                # maf_dict
                maf_dict[chr_num[:len(chr_num)-1]].setdefault(variant, minor_allele_frq)

    return maf_dict