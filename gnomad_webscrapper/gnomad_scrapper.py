
from selenium import webdriver
import pandas as pd
import argparse
import os
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By


class soup_parser():
    '''This class will take the html site and create a soup'''

    def __init__(self, variant_file_path: str, output_name: str):
        self.url_start = "https://gnomad.broadinstitute.org/variant/"

        self.url_end = "?dataset=gnomad_r2_1"

        self.var_list = self.get_var_list(variant_file_path)

        self.output_file_name = output_name

        self.delete_log()

        gnomad_dict: dict = self.get_freq()

        self.dict_to_csv(gnomad_dict)

        self.close_browser()

    def get_var_list(self, variant_file_path: str) -> list:
        '''This function will pull out the different rs values'''

        if variant_file_path[-4:] == ".csv":

            var_df: pd.DataFrame = pd.read_csv(variant_file_path, sep=",")

        # check if it is an excel file
        elif variant_file_path[-5:] == ".xlsx":

            var_df: pd.DataFrame = pd.read_excel(variant_file_path)
            print(var_df)
        # check if it is an excel file
        elif variant_file_path[-4:] == ".txt":

            var_df: pd.DataFrame = pd.read_csv(variant_file_path)

        var_list: list = var_df["RS Name"].values.tolist()
        print(len(var_list))
        return var_list

    @staticmethod
    def delete_log():
        if os.path.isfile('geckodriver.log'):

            os.remove('geckodriver.log')

    def get_freq(self) -> dict:

        browser = webdriver.Firefox()

        # creating a dictionary to keep track of the allele frequencies
        gnomad_freq_dict: dict = dict()

        # iterating through each variant
        for variant in self.var_list:

            # Getting the full page path for a specific variant
            full_url: str = "".join([self.url_start, variant, self.url_end])

            browser.get(full_url)

            # finding the table element
            print(browser.title)

            element = WebDriverWait(browser, 10).until(
                EC.presence_of_element_located(
                    (By.CLASS_NAME, "Table__BaseTable-sc-7fgtt2-0.PopulationsTable__Table-yt4zj1-0.gRZyOM"))
            )

            frequencies_table = browser.find_element_by_class_name(
                "Table__BaseTable-sc-7fgtt2-0.PopulationsTable__Table-yt4zj1-0.gRZyOM")

            pop_freq_txt: str = frequencies_table.text

            # This section gets the first index
            start_indx: int = pop_freq_txt.find("(non-Finnish)")

            start_indx = start_indx + len("(non-Finnish)")

            # find second index
            end_indx: int = pop_freq_txt.find("A", start_indx)

            allele_freq: str = pop_freq_txt[start_indx:end_indx].split(" ")[
                3].strip("\n")

            # adding the variant and the allele frequencies to a dictionary

            gnomad_freq_dict[variant] = allele_freq

        # closign the webdriver
        self.browser = browser

        return gnomad_freq_dict

    def close_browser(self):
        self.browser.quit()

    def dict_to_csv(self, dict: dict):

        # convert dictionary to dataframe
        gnomad_df = pd.DataFrame(dict.items(), columns=["RS Name", "MAF"])
        print(gnomad_df)
        gnomad_df.to_csv(self.output_file_name, index=False, sep="\t")


def run(args):
    "function to run"

    soup_parser(args.var_file, args.output)


def main():
    parser = argparse.ArgumentParser(
        description="")

    parser.add_argument("-v", help="This argument provides the file path to a file that list the variants, or a csv/excel file that has a specific column with the RS values",
                        dest="var_file", type=str, required=True)
    parser.add_argument("-o", help="This argument provides output path for the file",
                        dest="output", type=str, required=True)

    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
