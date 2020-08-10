#!/usr/bin/python

import sys  # THese are modules used
import argparse
import gzip
import multiprocessing as mp

from IBDinput_byBP_classes import Shared_Segment_Convert

####################################################################################################
# Initializing the parser


def run(args):
    print("running")
    for variant__info_tuple in zip(args.pheno, args.bp):
        pheno_file = variant__info_tuple[0]

        variant_position = variant__info_tuple[1]

        ibd_file_converter = Shared_Segment_Convert(
            args.input, pheno_file, args.output, args.format, args.min, args.thread, variant_position)

        parameter_dict = ibd_file_converter.generate_parameters()

        uniqID, dupID = ibd_file_converter.build_id_pairs()

        IBDdata, IBDindex = ibd_file_converter.create_ibd_arrays()

        ibd_file_converter.run_parallel(
            IBDdata, IBDindex, parameter_dict, uniqID)


def main():
    parser = argparse.ArgumentParser(
        description="This script converts the output of the  IBD detection programs to a human readable form")

    parser.add_argument("-i", '--input', help="This argument just list the path to the IBD detection software output",
                        dest="input", type=str, required=True)

    parser.add_argument("-p", "--pheno", help="This argument just list the path to a txt file of all the IIDs being focused on",
                        dest="pheno", nargs="+", type=str, required=True)

    parser.add_argument("-o", "--output", help="This argument list the output directory",
                        dest="output", type=str, required=True)

    parser.add_argument("-f", "--format", help="This argument specifies the IBD program used",
                        dest="format",  type=str, required=True)

    parser.add_argument("-m", "--min", help="This argument specifies the minimum cM threshold",
                        dest="min", type=str, required=True)

    parser.add_argument("-t", "--thread", help="This argument specifies the threadcount to be used in parallel",
                        dest="thread", type=str, required=True)

    parser.add_argument(
        "-b", "--bp", help="This argument specifies the nucleotide position of the variant of interest", dest="bp", nargs="+", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
