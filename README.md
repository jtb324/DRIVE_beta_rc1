# MEGA_CLI_project

This is a repository for working on turning the MEGA project into a commandline tool.

## Goal

The goal of this repository is to create a commandline tool that takes a raw file as input and then returns a list of all individual ids who contain a desired variant.

## Current Development

#TODO: Need to figure out how to make the multiVariantAnalysis run faster but trying to incorporate is in versus the nested for loops that use enumerate.
#TODO: Need to incorporate a function in the CLI to search through Pedigrees.

## Files in the Directory

- **MEGA_ID.py:** contains the argparse script to make the CLI. At the current development the program uses recodeA files from PLINK. Contains functions to determine the total # of individuals containing at least one variant and a function to find individuals containing multiple variants

- **IIDFindFunction.py:** contains five major functions to determine the number of individuals with at least one variant and the number of individuals with multiple variants.

  - totalVariantID: This function finds the total number of individuals carrying at least one variant. It prints the number of individuals carrying at least one variant to the console and then it creates a file "totalVariantIDList.txt" that contains the IIDs of each individual.

  * writePath: This simple function just creates a path to write files to. This was doing to keep the DRY principle of coding.

  * individualCount: This function forms the dictionary where the keys are the index of each variant in each row from the recode file. The values are the number of individuals carrying that combination of variants. This function then uses the csvDictWriter function to write this dictionary to a csv file.

  * multiVariantAnalysis: This function creates a dictionary where the keys are the index of each variant in each row from the recode file. The values are list of IIDs of individuals who carry that combination of variants. The csvDictWriter function is then used to write this dictionary to a csv file.

  * csvDictWriter: This function creates a csv file from dictionaries passed to it. It uses the write path function to create a path for the csv file to be written to. This function was also made to keep the DRY principle.

* **FindIndPedigree.py:** attempt at iterating through the list of total variants and seeing if any of those are in the Pedigree. Currently it is not writing properly. This is a working progress
