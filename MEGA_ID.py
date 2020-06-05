#! /usr/bin/env python

import argparse
from IdentifyingID import IdentifyID
from FindIndPedigree import findPedigreeInd


def run(args):
    if args.variantid == True:
        print("generating output....")
        print('\n')
        IdentifyID(args.input, args.output)

    elif args.pedigreeid == True:
        print("Generate the list of individuals in the pedigree files...")
        print('\n')
        findPedigreeInd(args.input, args.output)  # need to add output later


def main():
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK")

    parser.add_argument("-i", "--input", help="this is the pathway for the PLINK recoded input file.",
                        dest="input", type=str, required=True)

    parser.add_argument("-o", "--output", help="This is the directory that text files containing the ids of the individuals who have desired variants will be written to.",
                        dest="output", type=str, required=True)

    parser.add_argument("--variantID", help="This tag is used to indicate that the user wants to find the variats of interest in the raw files from PLINK",
                        dest="variantid", type=bool, default=False)

    parser.add_argument("--pedigreeID", help="This will find the matching values that are in the variant list formed by the variantID function from in the pedigree .genome_network ",
                        dest="pedigreeid", type=bool, default=False)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
