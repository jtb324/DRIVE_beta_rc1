# Import statements
import pandas as pd
import sys

######################################


class Pop_Filter:

    def __init__(self, pop_info_file: str, recoded_file):

        self.info_file = pop_info_file
        self.recode_file = recoded_file

    def load_files(self):
        '''This function loads both the file containing population info for the database, and the extracted recode file to pandas dataframes. The pop_info_file needs to contain a column titled "Pop" that has the population codes for each grid. This file also needs to have a column named ']"grid" with the IIDs for each grid.'''

        pop_info_df = pd.read_csv(self.info_file, sep="\t")
        print(self.recode_file)
        if isinstance(self.recode_file, pd.DataFrame):

            recode_df = self.recode_file
        else:

            recode_df = pd.read_csv(self.recode_file, sep=" ")
        print(recode_df)
        return pop_info_df, recode_df

    def get_pop_info_subset(self, pop_info_df, pop_code):
        '''This function subset the pop_info_df into a dataframe containing only specific population information'''

        try:
            pop_info_subset_df = pop_info_df[pop_info_df.Pop.isin([pop_code])]

        except KeyError:
            print("Make sure the the pop_info_file has a column named 'Pop'.")
            sys.exit(1)

        return pop_info_subset_df

    def filter_recode_df(self, pop_info_subset_df, recode_df):
        '''This function will filter the recoded file for only those where the IID is in the grid column.'''

        grid_list = pop_info_subset_df.grid.values.tolist()

        recode_df_filtered = recode_df[recode_df.IID.isin(grid_list)]

        return recode_df_filtered
