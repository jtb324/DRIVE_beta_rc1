#! /usr/bin/env python

import argparse
from IIDFindFunction import totalVariantID, multiVariantAnalysis, singleVariantAnalysis
from SearchPedigree import searchPedigree  # , multiCarriers


def run(args):
    if args.analysis == "total":
        print("generating list of IIDs who carry a desired variant....")
        totalVariantID(args.input, args.output)

    elif args.analysis == "multi":
        print("generating list of individuals carrying multiple variants....")
        multiVariantAnalysis(args.input, args.output, 'MultiVariantList.csv')

    elif args.analysis == "single":
        print("generating list of individuals at each variant index")
        singleVariantAnalysis(args.input, args.output, 'SingleVariantList.csv')

    elif args.analysis == "matchPED":
        print("generating a csv file of individuals found within the Pedigree...")
        searchPedigree(args.input, args.output,
                       args.drop_var, 'IndividInPedigree.csv')

    elif args.analysis == "multiIndivid":
        print("generating a csv file of individuals in the same pedigree who also are carrying the same variant")
        multiCarriers(args.input, args.output, 'multiIndivids.csv')


def main():
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK")

    parser.add_argument("--input", help="This is the pathway for the PLINK recoded input file. If you use the matchPED analysis argument then you should provide two input paths. The first is to the list of all variants. This should be a csv file and will have a list of the variant index and then a list of individuals who carry that variant. The second path is to the Pedigree file.At current development, this should be a .fam file",
                        dest="input", nargs="+", type=str, required=True)

    parser.add_argument("--output", help="This is the directory that text files containing the ids of the individuals who have desired variants will be written to.",
                        dest="output", type=str, required=True)

    parser.add_argument("--analysis", help="This tag indicates that the multiVariantAnalysis function will be called to analyze how many individuals carry multiple variants. Two csv files are made which contain the indices of the variants and a list of the individuals that contain those variants. This accepts single, total, multi, matchPED, multiIndivid", dest="analysis", type=str, default=False)

    parser.add_argument("--drop_var", help="This functionality is used to drop variants from a file if needed to for some reason. This is passed into the searchPedigree function incase maybe a certain variant is too common and can be removed",
                        dest="drop_var", type=str, nargs="+")

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
