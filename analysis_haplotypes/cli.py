import argparse
import analysis_haplotypes


def run(args):
    '''This is the main run function'''
    analysis_haplotypes.compare_haplotypes(args.hfile, args.threads,
                                           args.output)


def main():
    parser = argparse.ArgumentParser(
        description="This cli is used to run the full haplotype analysis")

    parser.add_argument("-o",
                        help="This argument list the path to output files in",
                        dest="output",
                        type=str,
                        required=True)

    parser.add_argument(
        "-hfile",
        help="This argument list the path to the haplotype_info__file",
        dest="hfile",
        type=str,
        required=True)

    parser.add_argument(
        "-t",
        help=
        "This argument list the number of threads to be used by the program",
        dest="threads",
        type=int,
        required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
