# MEGA_CLI_project

This is a repository for working on turning the MEGA project into a commandline tool.

## Goal

The goal of this repository is to create a commandline tool that takes a raw file as input and then returns a list of all individual ids who contain a variant.

## Current Development

- **MEGA_ID.py:** contains the argparse script to make the CLI. At the current development the program uses recodeA files from PLINK.

- **IdenifyingID.py:** contains the main function that opens the input file and reads it line by line. It then passes each line to the individualVariantID and the multiVariantID functions to identify individuals who have a single variant or those that have multiple variants. This file also passes an empty list to each of the above mentioned functions that the IDs of each variant carrying individual gets written to.

  - This file writes the ids of all the individuals with variants to a file named totalVariantIDList.txt. It then writes the ids of all the individuals with multiple variants to a file named multiVariantIDList.txt. The user provides the output directory and then these files are written into that directory.

  * Added a functionality where you can choose if you want to just find the ids in the raw file or use the raw file to search through provided Pedigrees. [EXPLAIN FURTHER]

- **VariantID.py:** contains the individualVariantID function that identifies individuals that contain a single variant. The IIDs for these individuals get appended to the totalVariantList that is passed as the second argument to the function.

- **MultipleVariantID.py:** contains the multiVariantID function that identifies individuals that contain multiple variants. These IIDs are also append to the multiVariantList that is passed as the second argument to the function.

* **FindIndPedigree.py:** attempt at iterating through the list of total variants and seeing if any of those are in the Pedigree. Currently it is not writing properly
