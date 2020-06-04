#! /usr/bin/env python

import argparse


def run(args):
    k = 5


def main():
    parser = argparse.ArgumentParser(
        description="This identifies individuals who have a specific variant in a raw file from PLINK")
    parser.add_argument("--input", help="this is the pathway for the PLINK recoded input file",
                        dest="input", type=str, required=True)
    parser.add_argument("--output", help="This is the pathway for the text file containing the ids of the individuals who have a desired variant",
                        dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
