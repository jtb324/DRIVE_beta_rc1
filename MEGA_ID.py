#! /usr/bin/env python

import argparse
from IIDFindFunction import totalVariantID, multiVariantAnalysis


def run(args):
    if args.totalvar == True:
        print("generating list of IIDs who carry a desired variant....")
        totalVariantID(args.input, args.output)

    if args.multivar == True:
        print("generating list of individuals carrying multiple variants....")
        multiVariantAnalysis(args.input, args.output, 'MultiVariantList.csv')

    # elif args.pedigreeid == True:
    #     print("Generate the list of individuals in the pedigree files...")
    #     print('\n')
    #     findPedigreeInd(args.input, args.output)  # need to add output later


def main():
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK")

    parser.add_argument("--input", help="this is the pathway for the PLINK recoded input file.",
                        dest="input", type=str, required=True)

    parser.add_argument("--output", help="This is the directory that text files containing the ids of the individuals who have desired variants will be written to.",
                        dest="output", type=str, required=True)

    parser.add_argument("--totalVar", help="this tag is used to find the total number of individuals containing a variant and outputs a txt file containing all the IIDS.",
                        dest="totalvar", type=bool, default=False)

    parser.add_argument("--multiVar", help="This tag indicates that the multiVariantAnalysis function will be called to analyze how many individuals carry multiple variants. Two csv files are made which contain the indices of the variants and a list of the individuals that contain those variants.", dest="multivar", type=bool, default=False)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
