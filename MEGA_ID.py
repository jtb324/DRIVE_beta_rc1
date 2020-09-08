#! /usr/bin/env python
import logging
import argparse
import os.path
from os import path
import pandas as pd

from identify_single_var_carrier import singleVariantAnalysis
from SearchPedigree import searchPedigree
from allele_count import allele_counts
from identify_multi_var_carriers import multiVariantAnalysis
from create_networks import create_networks
from csv_writer_class import Csv_Writer_Object


def run(args):

    if args.analysis == "multi":

        log_format = '%(asctime)s - %(levelname)s : %(message)s'

        logging.basicConfig(
            filename="".join([args.output, '/multi_variant_analysis.log']), level=logging.INFO,
            format=log_format)

        logging.info(
            'Generating a list of IIDs who carry multiple variants...')

        print("generating list of individuals carrying multiple variants....")

        multiVariantAnalysis(args.input, args.output,
                             args.compatible_format, 'multi_variant_list.csv')

        logging.info(
            'Finished creating a list of individuals carrying multiple variants')

    elif args.analysis == "single":

        log_format = '%(asctime)s - %(levelname)s : %(message)s'

        logging.basicConfig(filename="".join([args.output,
                                              '/single_variant_analysis.log']), level=logging.INFO,
                            format=log_format)

        logging.info("generating list of individuals at each probe id...")

        print("generating list of individuals at each probe id...")

        # The args.input should be a directory indicating where the raw files are located
        singleVariantAnalysis(args.input, args.output, args.compatible_format,
                              'single_variant_list.csv', args.pop_info, args.pop_code)

        logging.info(
            'Finished creating a file of individuals for each probe id')

    elif args.analysis == "matchPED":

        log_format = '%(asctime)s - %(levelname)s : %(message)s'

        logging.basicConfig(filename="".join([args.output,
                                              '/search_pedigree_analysis.log']), level=logging.INFO,
                            format=log_format)

        logging.info(
            "identifying carriers within the provided pedigree,{}". format(args.output[1]))

        print("generating a csv file of individuals found within the Pedigree...")

        searchPedigree(args.input, args.output,
                       args.drop_var, args.compatible_format, args.pedigreeSubset, 'all_ind_in_pedigree.csv')

    elif args.analysis == "allele_counts":

        log_format = '%(asctime)s - %(levelname)s : %(message)s'

        logging.basicConfig(filename="".join([args.output,
                                              '/allele_count_analysis.log']), level=logging.INFO,
                            format=log_format)

        logging.info(
            "determining the allele_counts each variant for families containing carriers within the provided network file...")

        print("generating list of the allele counts for each network...")

        allele_counts(args.input, args.fam_file, args.output)

    elif args.analysis == "draw_networks":

        print("generating pdf files of networks of individuals who share segments...")

        # This dictionaru keeps track of how many carriers are actually in the network. It needs to be a global variable so that it is just extended for each variant instead of recreated
        carrier_in_network_dict = dict()

        # being able to draw networks for numerous variants
        variant_info_df = pd.read_csv(args.var[0], sep=" ", header=None, names=[
            "output_file_name", "variant_bp", "variant_id"])

    # TODO: get rid of this loop that gathers everything. This will instead be done in the actually create_networks function
    # no longer need to get the different value for a txt file. can extract the variant_ids from teh allpair.new.txt files

    # can still pass hte output. But the function will generate the variant name

        for info_tuple in zip(args.segments_file, file_name_list, variant_id_list):

            print(info_tuple)
            segment_file = info_tuple[0]
            output_file_name = info_tuple[1]
            var_of_interest = info_tuple[2]

            variant_file = args.var_file

            output = "".join([args.output, "network_imgs"])

            if not path.exists(output):
                os.mkdir(output)

            log_format = '%(asctime)s - %(levelname)s : %(message)s'

            logging.basicConfig(filename="".join([output,
                                                  '/drawing_networks.log']), level=logging.INFO,
                                format=log_format)

            logging.info(
                "Drawing networks for all individuals identified as shared segments...")

            carrier_in_network_dict = create_networks(segment_file, variant_file, carrier_in_network_dict,
                                                      var_of_interest, output)

        logger = logging.getLogger(args.output+'/carriers_in_network.log')

        # Writing the dictionary to a csv file
        csv_writer = Csv_Writer_Object(
            carrier_in_network_dict, args.output, "carriers_in_networks.csv", logger)

        csv_writer.log_file_path()

        csv_writer.write_to_csv()


def main():
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK")

    parser.add_argument("--input", help="This is the pathway for the PLINK recoded input file. If you use the matchPED analysis argument then you should provide two input paths. The first is to the list of all variants. This should be a csv file and will have a list of the variant index and then a list of individuals who carry that variant. The second path is to the Pedigree file.At current development, this should be a .fam file.",
                        dest="input", nargs="+", type=str)

    parser.add_argument("--output", help="This is the directory that text files containing the ids of the individuals who have desired variants will be written to.",
                        dest="output", type=str, required=True)

    parser.add_argument("--analysis", help="This tag indicates that the multiVariantAnalysis function will be called to analyze how many individuals carry multiple variants. Two csv files are made which contain the indices of the variants and a list of the individuals that contain those variants. This accepts single, total, multi, matchPED, allele_counts, draw_networks", dest="analysis", type=str, default=False)

    parser.add_argument("--drop_var", help="This functionality is used to drop variants from a file if needed to for some reason. This is passed into the searchPedigree function incase maybe a certain variant is too common and can be removed",
                        dest="drop_var", type=str, nargs="+")

    parser.add_argument("--format", help="This argument will enable several additional output files to be created that are easier for non python programs to interact with",
                        dest="compatible_format", type=bool, default=False)

    parser.add_argument("--pedigreeSubset", help="This option can be set so that the output either returns a list of individuals filtered by if the IID does not equal the FID or if it returns every individual found within the provided network file. By selecting '--pedigreeSubset Full' then it will list the full size of the pedigree.", dest="pedigreeSubset", type=str, default=False)

    parser.add_argument("--fam_file", help="This provides a path to the desired .fam of networks. This argument is used in the allele_count script",
                        dest="fam_file", type=str, default=False)

    parser.add_argument("--shared_segments_file", help="This argument provides a path to the list of shared segments that can be used to form networks",
                        dest="segments_file", type=str, nargs="+", default=False)

    parser.add_argument("--variant_file", help="This argument provides a path to a file that list all individuals that carry a specific variant",
                        dest="var_file", type=str, default=False)

    parser.add_argument("--var_of_interest", help="This argument passes a variant of interest that can filter down dataframes when trying to draw networks.",
                        dest="var", type=str, nargs="+", default=False)

    parser.add_argument("--pop_info", help="This argument provides the file path to a file containing the population distribution of a dataset for each grid. This file should be a text file and at least contain two columns, where one column is 'Pop', the population code for each grid based on 1000 genomes, and then the second column is 'grid', which list the IIDs for each each grid.", dest="pop_info", type=str, required=False)

    parser.add_argument("--pop_code", help="This argument can be a single population code or a list of population codes that someone is looking for. Population codes have to match the 1000 genomes.",
                        dest="pop_code", type=str, required=False, nargs="+")

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
