# MEGA_CLI_project

This is a repository for working on turning the MEGA project into a commandline tool using the argparse package in python.

## Goal

The goal of this repository is to create a commandline tool that takes a raw file as input and then returns a list of all individual ids who contain a desired variant.

## Current Development

## Files in the Directory

- **shell_script**: This is a folder containing shell script I run. The ∏keyidea behind this was it may help with reproducibility because I can put what files I used for input, what the output path was, etc..

  - **MEGA_ID_sh_script.sh**: This was a script used to run the matchPED analysis in the MEGA_ID.py program. It took two file inputs. The first input was file path to the csv file containing two columns, where one is MEGA_ex array variant id and then the other column is a list of IIDs who carry that variant. The second input is a pedigree file. This was run using a .fam file. The the output directory is specified. The analysis type was "matchPED" which is used to identify all individuals in the provided pedigree file and pairs them to the corresponding variant id. This analysis type also provides other information about which indiivduals are in the network and then it determines if there are multiple individuals in the same network who also carry the same variant.

* **MEGA_ID.py:** contains the argparse script to make the CLI. At the current development the program uses recodeA files from PLINK. Contains functions to determine the total # of individuals containing at least one variant and a function to find individuals containing multiple variants

* **IIDFindFunction.py:** contains five major functions to determine the number of individuals with at least one variant and the number of individuals with multiple variants.

  - totalVariantID: This function finds the total number of individuals carrying at least one variant. It prints the number of individuals carrying at least one variant to the console and then it creates a file "totalVariantIDList.txt" that contains the IIDs of each individual.

  * writePath: This simple function just creates a path to write files to. This was doing to keep the DRY principle of coding.

  * individualCount: This function forms the dictionary where the keys are the index of each variant in each row from the recode file. The values are the number of individuals carrying that combination of variants. This function then uses the csvDictWriter function to write this dictionary to a csv file.

  * multiVariantAnalysis: This function creates a dictionary where the keys are the index of each variant in each row from the recode file. The values are list of IIDs of individuals who carry that combination of variants. The csvDictWriter function is then used to write this dictionary to a csv file.

  * csvDictWriter: This function creates a csv file from dictionaries passed to it. It uses the write path function to create a path for the csv file to be written to. This function was also made to keep the DRY principle.

- **SearchPedigree.py**: This script contains a function to search through a provided .fam pedigree file. The function takes both a list of Variant IIDs and the .fam pedigre family as inputs as well as a filename for the output file.

  There are several functions in it:

  - writePath: same as the IIDFindFunction.py

  - csvDictWriter: the same as IIDFindFunction.py

  * drop_variant: This function can allow the used to drop specified variants from the list of input Variant IIDs. This was added incase the prior functions in the IIDFindFunction.py and indicate that there are variants that are more common than expected.

  - pedigreeCount: This function searches through the provided dictionary and will determine the number of individuals carrying each variant.

  * network_sizes: This function returns some information about the network size. It first returns information about the number of individuals within each network found within the pedigree. This does not necessarily indicate the size of the network. It just list the number of individuals found per network. This step returns a file called network_counts.csv. The next step just list all the individuals found per network. It returns a csv file, titled network_list.csv, where one column is the FID which is the network file and the other column is a list of IIDs for each carrier in the variant.

  - searchPedigree: This function matches IID from the provided variant file to IIDs in the provided pedigree and creates a csv file containing a list of IIDS for each variant, if the IID does not equal the FID.

  * multi_ind_in_pedigree: This function accepts three arguments. The first input is a dictionary formed in the searchPedigree function, where the key is the variant id and the value is a list of individual carriers. This dictionary basically list all individuals found within the .fam for a specific variant, but it does not tell anything about the networks. The second one is a subset of the original pedigree dataframe in the searchPedigree function that contains only IIDs that were found in the variant list and in the pedigree. The third argument is a path to write the output file to. In this function the lists of a individuals for a specific variant, from the first argument, are used, individually, to further subset the pedigree subset so that the dataframe on contains those specific individuals in the list. This subset is then grouped by the network "FID" and they are counted so that you can see how many individuals are in each small subset. All networks where there is only one individual are then filtered out. The FIDs of the remaining networks are writing to a list are used to specifically find all the IIDs in the original pedigree dataframe from the second argument. It then returns a two column csv file, titled multiple_ind_in_pedigree.csv, where the first column is a tuple containing the variant id and then the network FID. The second column contains a list individuals who carry the same variant in the same network.

* **NetworkSize.py**: This is an attempt to determine the size of each network that was identified in the searchPedigree function. The current attempt takes the FID list from the ind_network_count.csv function and opens the file in the correct directory and then will try to find all the unique individuals in the network file. Need to try to read in file in pandas.
