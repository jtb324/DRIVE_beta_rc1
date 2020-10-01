
from selenium import webdriver
import pandas as pd
import argparse
import os
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from webdriver import determine_webdriver
from selenium.webdriver.firefox.options import Options
import sys
import re


class soup_parser():
    '''This class will take the html site and create a soup and 
    then extract hte deired minor allele frequency'''

    def __init__(self, variant_file_path: str, output_name: str, browser_name: str, wait_time: int):
        self.url_start: str = "https://gnomad.broadinstitute.org/variant/"

        self.url_end: str = "?dataset=gnomad_r2_1"

        self.output_file_name: str = output_name

        self.browser: str = browser_name.lower()

        self.time: int = wait_time

        self.var_list: list = self.get_var_list(variant_file_path)

        self.delete_log()

        gnomad_dict: dict = self.get_freq()

        self.dict_to_csv(gnomad_dict)

    def get_var_list(self, variant_file_path: str) -> list:
        '''This function will pull out the different rs values'''

        if variant_file_path[-4:] == ".csv":

            var_df: pd.DataFrame = pd.read_csv(variant_file_path, sep=",")

        # check if it is an excel file
        elif variant_file_path[-5:] == ".xlsx":

            var_df: pd.DataFrame = pd.read_excel(variant_file_path)

        # check if it is an excel file
        elif variant_file_path[-4:] == ".txt":

            var_df: pd.DataFrame = pd.read_csv(variant_file_path)

        var_list: list = var_df["RS Name"].values.tolist()

        return var_list

    @staticmethod
    def delete_log():
        '''This method removes the geckodriver.log if it is present from
        a previous run'''

        if os.path.isfile('geckodriver.log'):

            os.remove('geckodriver.log')

    def get_freq(self) -> dict:
        '''This function gets the minor allele frequencies for variants of interest'''

        # creating a headless mode
        options = Options()

        options.headless = True
        # Creating a dictionary to act as a handler based on what browser is used

        browser_handling_dict: dict = {
            "firefox": determine_webdriver(self.browser, options),
            "chrome": determine_webdriver(self.browser, options),
        }

        # creating the driver object
        browser = browser_handling_dict[self.browser]

        # creating a dictionary to keep track of the allele frequencies
        gnomad_freq_dict: dict = dict()

        # iterating through each variant
        for variant in self.var_list:

            # Getting the full page path for a specific variant
            full_url: str = "".join([self.url_start, variant, self.url_end])

            browser.get(full_url)

            # finding the table element
            element = WebDriverWait(browser, self.time).until(
                EC.presence_of_element_located(
                    (By.CLASS_NAME, "Table__BaseTable-sc-7fgtt2-0.PopulationsTable__Table-yt4zj1-0.gRZyOM"))
            )

            frequencies_table = browser.find_element_by_class_name(
                "Table__BaseTable-sc-7fgtt2-0.PopulationsTable__Table-yt4zj1-0.gRZyOM")

            # Getting the text from the table element
            pop_freq_txt: str = frequencies_table.text

            # This section gets the first index
            start_indx: int = pop_freq_txt.find("(non-Finnish)")

            start_indx = start_indx + len("(non-Finnish)")

            # Get a substring that starts at the end of the non-Finnish
            pop_freq_substring: str = pop_freq_txt[start_indx:]

            # Use a regular expression to find The first letter in the substring
            letter_match = re.search("[a-zA-Z]", pop_freq_substring)

            first_letter: str = letter_match.group(0)

            # find second index
            end_indx: int = pop_freq_substring.find(first_letter)

            allele_freq: str = pop_freq_substring[:end_indx].split(" ")[
                3].strip("\n").strip(" ")

            print(
                f"This is the allele frequencies for the variant {variant} is {allele_freq}")

            # adding the variant and the allele frequencies to a dictionary

            gnomad_freq_dict[variant] = allele_freq

        # closign the webdriver
        browser.close()
        browser.quit()

        return gnomad_freq_dict

    def dict_to_csv(self, dict: dict):

        # convert dictionary to dataframe
        gnomad_df = pd.DataFrame(dict.items(), columns=["RS Name", "MAF"])

        # writing the dataframe to a file
        gnomad_df.to_csv(self.output_file_name, index=False, sep="\t")


def run(args):
    "function to run"
    # Checks to make sure that the browser being used is a supported browser
    if args.browser.lower() not in ["chrome", "firefox"]:
        print("please use one of the two supported web browser: Firefox or Google Chrome")
        sys.exit(1)

    soup_parser(args.var_file, args.output, args.browser, args.wait_time)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("-v", help="This argument provides the file path to a file that list the variants, or a csv/excel file that has a specific column with the RS values",
                        dest="var_file", type=str, required=True)

    parser.add_argument("-o", help="This argument provides output path for the file",
                        dest="output", type=str, required=True)

    parser.add_argument("-b", help="This argument provides name of the internet browser used",
                        dest="browser", type=str, required=True)

    parser.add_argument("-t", help="This parameter gives a amount of time in seconds that the webdriver should wait for a http response. This defaults to 20",
                        dest="wait_time", default=20, type=int, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
