## [beta 1.0.0 rc1] - 2021-06-2:

### ##[Unreleased]

### ##Added:

### ##Changed:

***What was done:***

- Updated lines 83 to 123 in the identify_single_var_carriers.py file
    - These updates included changing the output to be a dictionary of dictionaries where the outer key is the chromosome number of the raw file being used, the inner key is the variant probe id, and the inner value is a list of IIDs for that variant. If there are no genotyped carriers for the variant then the list will contain a single value of N/A.

### ##TODO:

***What needs to be done:***

- Need to create a function that can write the dictionary to an output file. should probably change these files from a csv file to a text file.
    - Consider making this into a single file that has a third column for the chromosome number
- Also need to return the dictionary for use in the next function

---

## [beta 1.0.0 rc1] - 2021-06-3:

### ##[Unreleased]

### ##Added:

***What was done:***

- added a function called write_to_file in the identify_single_var_carriers.py file. This function iterates through the dictionary to pull out the iids that carry the variant of interest for each variant for each chromosome.
- Added a change to the get_allele_freq function that returns a dictionary that has the minor allele frequency for each variant for each chromosome. This dictionary is the maf_dict

### ##Changed:

- switch the output file from a csv file to a tab delimited text file
- Add a try and except statement to the get_chr_num in the get_chr_num.py file.
    - This except catches the Attribute error that would be raised if there was no match for the regular expression
        - If there is no match then the get_chr_num function will return a string value of "0"
        - If there is a match then the function returns a chromosome number of the format "chrXX" where XX is a number
- Added a progress bar to the single_variant_analysis function using the tqdm package
    - This just visually helps the user understand what is happening.
    - Added this dependency to the requirements.txt file
- refactored the get_allele_frq function in the maf_determination.py file.
    - Accepts the carrier_dict as an input. This dictionary is the returned result of the single_variant_analysis function from the identify_single_var_carriers.py which contains the list of carriers per variant per chromosome.
    - Also made the function accept the string of the directory that has the raw files instead of the list of the raw files. The function will now create the list within the function
    - Adjusted the function so that it iterates through the variants for each chromosome and uses the value_counts() from the pandas series to count the number of heterozygous and homozygous carriers.
        - These value counts are used to determine the allele count for the minor allele which is then used to determine the minor allele frequency.

### ##Removed:

- removed all the sections that combined the different files into a single dataframe per chromosome
    - This structure was replaced by a dictionary structure mentioned in Jun 2, 2021 changelog.
    - This dictionary is written to a single file in the new added function
    - This dictionary contains information about the carriers for each variant for every chromosome. This approach is different then the old approach where individual files were created for each chromosome
- Removed the code that created a different csv file for each chromosome. Instead they should all be in one text file now
- Removed the determine_maf function because it was redundant.
    - Instead the program now just uses the get_allele_frq function

### ##TODO:

- Need to fix the unit test for the determine_maf.py script
- Need to move on the check_maf function to refactor it to use the maf_dict to check for minor allele frequencies thresholds