# MEGA_CLI_project

This is a repository for working on turning the MEGA project into a commandline tool using the argparse package in python.

## Goal

The goal of this repository is to create a commandline tool that takes a raw file as input and then returns a list of all individual ids who contain a desired variant.

## Current Development

#TODO: Need to figure out how to make the multiVariantAnalysis run faster but trying to incorporate is in versus the nested for loops that use enumerate.
#TODO: Need to determine a way to identify if there are multiple individuals in the same pedigree who carry the same variants.

## Files in the Directory

- **MEGA_ID.py:** contains the argparse script to make the CLI. At the current development the program uses recodeA files from PLINK. Contains functions to determine the total # of individuals containing at least one variant and a function to find individuals containing multiple variants

- **IIDFindFunction.py:** contains five major functions to determine the number of individuals with at least one variant and the number of individuals with multiple variants.

  - totalVariantID: This function finds the total number of individuals carrying at least one variant. It prints the number of individuals carrying at least one variant to the console and then it creates a file "totalVariantIDList.txt" that contains the IIDs of each individual.

  * writePath: This simple function just creates a path to write files to. This was doing to keep the DRY principle of coding.

  * individualCount: This function forms the dictionary where the keys are the index of each variant in each row from the recode file. The values are the number of individuals carrying that combination of variants. This function then uses the csvDictWriter function to write this dictionary to a csv file.

  * multiVariantAnalysis: This function creates a dictionary where the keys are the index of each variant in each row from the recode file. The values are list of IIDs of individuals who carry that combination of variants. The csvDictWriter function is then used to write this dictionary to a csv file.

  * csvDictWriter: This function creates a csv file from dictionaries passed to it. It uses the write path function to create a path for the csv file to be written to. This function was also made to keep the DRY principle.

* **SearchPedigree.py**: This script contains a function to search through a provided .fam pedigree file. The function takes both a list of Variant IIDs and the .fam pedigre family as inputs as well as a filename for the output file.

There are several functions in it:

- writePath: same as the IIDFindFunction.py

- csvDictWriter: the same as IIDFindFunction.py

- pedigreeCount: This function searches through the provided dictionary and will determine the number of individuals carrying each variant.

- searchPedigree: This function matches IID from the provided variant file to IIDs in the provided pedigree and creates a csv file containing a list of IIDS for each variant, if the IID != equal the FID.

* multiCarriers: this function determines if there are multiple individuals who carry one variant in the same pedigree. It outputs a csv file. This function is still in development and requires the use of "multiIndivid" as the argument in the analysis flag.
