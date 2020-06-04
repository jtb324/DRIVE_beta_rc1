# MEGA_CLI_project

This is a repository for working on turning the MEGA project into a commandline tool

## Goal

The goal of this repository is to create a commandline tool that takes a raw file as input and then returns a list of all individual ids who contain a variant

## Current Development

- MEGA_ID.py: contains the argparse script to make the CLI
- IdenifyingID.py: contains the main function that opens the input file and reads it line by line. It then passes each line to the individualVariantID and the multiVariantID functions to identify individuals who have a single variant or those that have multiple variants.
- VariantID.py: contains the individualVariantID function that identifies individuals that contain a single variant
- MultipleVariantID.py: contains the multiVariantID function that identifies individuals that contain multiple variants.
