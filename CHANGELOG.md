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
- Added logging to the load_pheno_file so that if the file is not found then it will log a message to the log file

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
* Changed the check_maf_range.py file. Now it only has two functions: check_maf and check_frequencies. check_maf function is the main one used by the program
    - Now the function utilizes the maf_dict to get the allele frequencies for each variant and check them against the threshold. This way the program does not have to open up unnecessary files.
    - Also was able to just take advantage of the filter function to filter the variants instead of using the three other functions.
    - The function still returns a tuple where the first value is a list of variants above the threshold and the second value is the "program_end_code" where a 1 means the user wishes to terminate the program and a 0 means the user wishes to continue with the program or none of the variants exceed the threshold.

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
___

## [beta 1.0.0 rc1] - 2021-06-6:

### ##[Unreleased]

### ##Added:

***What was done:***

- Created a dataclass object called Readme_info in the file_generator.py file.
    - This dataclass object has two attributes called the readme_output_path and the read_body_text where the first one is a string with the output path to write the readme to and the second attribute is either a string or a list of strings containing the text that will be written to the README

### ##Changed:

### ##Removed:

- Removed the carrier_dir parameter in both the collect_files function from the gather_ibd_info.py file and from the convert_ibd_func_param from the DRIVE.py file
- Removed the carrier_file_list input value from the get_file_dict function in the create_file_dict.py function.
- also removed the carrier key from the find_match function in the create_file_dict.py function since it was removed in the above to points
- removed the readme_class.py file and the other file from the readme_generators module but kept the readme_body_text.py.
    - This was done because the code is repeated in the file_generator.py files
### ##TODO:
___
## [beta 1.0.0 rc1] - 2021-06-9:

### ##[Unreleased]

### ##Added:

***What was done:***

- Added some functions to the Readme_Info that will automatically create the readme when the object is created.
- Added two modules: This might be restructure to reflect the gene or phenotype driven approaches
    - **One was carrier_identification**
    - **One was segment analysis**
    - These two modules contain all the scripts that are dependent on these parts.

### ##Changed:

- decided to begin restructuring the project so that the cli can be run from three commands. This is being done to simplify the [DRIVE.py](http://drive.py) file.
    - The first command will determine which individuals carry the variants of interest.
        - It will need to also potentially run the plink step and determine the minor allele frequencies
    - The second command will do all the shared segment analysis and will create the networks
    - The third command will get all the information for each individual if the user wants that
- Changed the add_line function in the Readme class in the file_generator.py file so that it can accept a string or a list. This was done incase the readme needs to have more than one line of info

### ##Removed:

### ##TODO:

- Need to fix the readmes for each section
___
## [beta 1.0.0 rc1] - 2021-06-11:

### ##[Unreleased]

### ##Added:

***What was done:***

- Finished the main function in the carrier_id_main.py file that is responsible for running the single_variant_analysis function.
- Added a function called plink_runner that handles all the calls to plink depending on what values the user passes to the initial determine carriers function
- Added option arguments such as:
    - run_plink
    - variant_file
    - binary_file
    - RANGE_STR
    - RANGE_END
- The above arguments are need to run plink if the user chooses to do so in the program. If the user does not need to do so then the values are None
- Also added a MAF_THRESHOLD value to the determine_carriers function as well as a RECODE_FLAGS value. These are also need to run plink and are only used in the if statement if the run_plink value is True (which it is by default)
- Added a dictionary to keep track of user parameters and write them to the log file

### ##Changed:

- Removed alot of the initial lines in the drive script that dealt with determine parameters and moved them into the carrier_id_main.py file as arguments for the determine_carriers cli function
- Also moved all of the functions that check the minor allele frequency to the carrier_id_main.py file in the main function

### ##Removed:

- Removed the need to have the readme filepath and the text as parameters to the single_variant_analysis function because the readme is just created when it is initialized in line 34 now.
- Change the cli in line 26 of the [DRIVE.py](http://drive.py) file to add the second cli from the carrier_id_main.py file as an option argument. This will also be done for the shared_segment detection scripts

### ##TODO:

- Need to fix the readmes for each section
___
## [beta 1.0.0 rc1] - 2021-06-14:

### ##[Unreleased]

### ##Added:

***What was done:***

- Added the segment_analysis_app into the main cli in the [DRIVE.py](http://drive.py) file. This way the user can select this analysis type as an option
    - It has the name segment_analysis in the user choices
- Added The main command function main into the segment_analysis_main.py file. This will have all the arguments that the user can enter such as:
    - output
    - pheno_gmap
    - pheno_carriers
    - MIN_CM
    - THREADS
    - MAF_THRESHOLD
    - IBD_programs

### ##Changed:

- moved the collect_phenotype_info module into the segment_analysis module  instead of in the carrier_identification module because it is used in the segment_analysis module.
- Changed the program so there will be a log file for the carrier_detection part and a log file for the segment_analysis_main.py part

### ##Removed:

- Started to remove the individual lines in the [DRIVE.py](http://drive.py) file that are being replaced in the segment_analysis_main.py file

### ##TODO:

- Need to fix the readmes for each section
- Need to branch the mainn function in the segment_analysis_main.py file so that it does one analysis for the phenotype side and one for the gene based analysis

___
## [beta 1.0.0 rc1] - 2021-06-16:

### ##[Unreleased]

### ##Added:

***What was done:***

- Added back in the file_dict_creator  because that script was actually used through the program

### ##Changed:

- 

### ##Removed:

- Started to remove the individual lines in the [DRIVE.py](http://drive.py) file that are being replaced in the segment_analysis_main.py file

### ##TODO:

- The program is no longer identifying the test set. This needs to be fixed

---

## [beta 1.0.0 rc1] - 2021-06-17:

### ##[Unreleased]

### ##Added:

***What was done:***
* Added a check to maKe sure that the inputs were correct and that the program checks the columns for the carrier file for the phenotype analysis
### ##Changed:

- Fixed the phenotype based analysis

* check all the logic for the phenotype based analysis and compared its output to the dcm output
    * The expected five carriers were identified in both the gene based analysis and the phenotype based analysis


### ##Removed:

### ##TODO:
* fix readmes for the formatted_ibd_output and the networks directory
* Work on CLI design so think about how to make like a progress bar to let people know it is still going
* work on making sure that the script can run in parallel
* work on reimplementing the part that gets the segment lengths
* fix unit test
    * at the moment every unit test should fail

## [beta 1.0.0 rc1] - 2021-06-23:

### ##[Unreleased]

### ##Added:

***What was done:***

### ##Changed:

- Fixed a bug where one one variant was working per chromosome for the gene based analysis. The dictionary was being created wrong so it was overwritting the dictionary every time



### ##Removed:

### ##TODO:
* fix readmes for the formatted_ibd_output and the networks directory
* Work on CLI design so think about how to make like a progress bar to let people know it is still going
* work on making sure that the script can run in parallel
* work on reimplementing the part that gets the segment lengths
* fix unit test
    * at the moment every unit test should fail
* Need to understand why for the phenotype null set that there are pairs where both iids are the same