import argparse
import analysis_haplotypes


def run(args):
    '''This is the main run function'''
    analysis_haplotypes.gather_haplotypes(args.hfile, args.output, args.bfile,
                                          args.pop_file, args.pop_code)


def main():
    parser = argparse.ArgumentParser(
        description="This cli is used to run the full haplotype analysis")

    parser.add_argument("-o",
                        help="This argument list the path to output files in",
                        dest="output",
                        type=str,
                        required=True)

    parser.add_argument(
        "-bfile",
        help="This argument list the path to the binary file used by plink",
        dest="bfile",
        type=str,
        required=True)

    parser.add_argument(
        "-hfile",
        help="This argument list the path to the haplotype_info__file",
        dest="hfile",
        type=str,
        required=True)

    parser.add_argument(
        "-pop_info",
        help="This argument list the file path to the population info file",
        dest="pop_file",
        type=str,
        required=True)

    parser.add_argument(
        "--pop_code",
        help="This argument list the population code to be used",
        dest="pop_code",
        type=str,
        required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
