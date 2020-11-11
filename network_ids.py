# this script will add the network ids to the file
import pandas as pd


def get_network_id(network_df_subset: pd.DataFrame) -> str:
    '''This function will get the network id for the subsetted dataframe'''

    network_id: str = network_df_subset["Network ID"].tolist()[0]

    return network_id


def filter_for_pairs(network_df: pd.DataFrame, pair_1: str, pair_2: str) -> set:
    '''This function will filter down the dataframe to just the ideally two rows that 
    contain information about the pairs'''

    network_df_subset: pd.DataFrame = network_df[(
        network_df.IID == pair_1) | (network_df.IID == pair_2)]

    network_id: str = get_network_id(network_df_subset)

    return network_id

# get the network id dataframe into the a panda's dataframe


def filter_df(network_file: str, variant_id: str) -> pd.DataFrame:
    '''This function will filter the dataframe to only the only that contains the specific variant'''

    network_df: pd.DataFrame = pd.read_csv(network_file)

    network_df_subset: pd.DataFrame = network_df[network_df.variant_id == variant_id]

    return network_df_subset
