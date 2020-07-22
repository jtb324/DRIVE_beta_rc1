#! /usr/bin/env python

import argparse
from IIDFindFunction import totalVariantID, multiVariantAnalysis, singleVariantAnalysis
from SearchPedigree import searchPedigree
from allele_count import allele_counts


def run(args):
    if args.analysis == "total":
        print("generating list of IIDs who carry a desired variant....")
        totalVariantID(args.input, args.output)

    elif args.analysis == "multi":
        print("generating list of individuals carrying multiple variants....")
        multiVariantAnalysis(args.input, args.output,
                             args.compatible_format, 'multi_variant_list.csv')

    elif args.analysis == "single":
        print("generating list of individuals at each variant index")
        singleVariantAnalysis(args.input, args.output, args.compatible_format,
                              'single_variant_list.csv')

    elif args.analysis == "matchPED":
        print("generating a csv file of individuals found within the Pedigree...")
        searchPedigree(args.input, args.output,
                       args.drop_var, args.compatible_format, args.pedigreeSubset, 'all_ind_in_pedigree.csv')

    elif args.analysis == "allele_counts":
        print("generating list of the allele counts for each network...")
        allele_counts(args.input, args.fam_file, args.output)


def main():
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK")

    parser.add_argument("--input", help="This is the pathway for the PLINK recoded input file. If you use the matchPED analysis argument then you should provide two input paths. The first is to the list of all variants. This should be a csv file and will have a list of the variant index and then a list of individuals who carry that variant. The second path is to the Pedigree file.At current development, this should be a .fam file.",
                        dest="input", nargs="+", type=str, required=True)

    parser.add_argument("--output", help="This is the directory that text files containing the ids of the individuals who have desired variants will be written to.",
                        dest="output", type=str, required=True)

    parser.add_argument("--analysis", help="This tag indicates that the multiVariantAnalysis function will be called to analyze how many individuals carry multiple variants. Two csv files are made which contain the indices of the variants and a list of the individuals that contain those variants. This accepts single, total, multi, matchPED, allele_counts", dest="analysis", type=str, default=False)

    parser.add_argument("--drop_var", help="This functionality is used to drop variants from a file if needed to for some reason. This is passed into the searchPedigree function incase maybe a certain variant is too common and can be removed",
                        dest="drop_var", type=str, nargs="+")

    parser.add_argument("--format", help="This argument will enable several additional output files to be created that are easier for non python programs to interact with",
                        dest="compatible_format", type=bool, default=False)

    parser.add_argument("--pedigreeSubset", help="This option can be set so that the output either returns a list of individuals filtered by if the IID does not equal the FID or if it returns every individual found within the provided network file. By selecting '--pedigreeSubset Full' then it will list the full size of the pedigree.", dest="pedigreeSubset", type=str, default=False)

    parser.add_argument("--fam_file", help="This provides a path to the desired .fam of networks. This argument is used in the allele_count script",
                        dest="fam_file", type=str, default=False)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
