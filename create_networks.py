import logging

from network_creator_class import Network_Img_Maker


def create_networks(segments_file, variant_file, ind_in_network_dict, variant_of_interest, output_path):

    logger = logging.getLogger(output_path+'/drawing_networks.log')

    network_drawer = Network_Img_Maker(
        segments_file, variant_file, output_path, logger)

    iid_list = network_drawer.isolate_variant_list(variant_of_interest)

    # this line loads the provided dataframe using a tab separator with no header value. it also skips the first row
    # which only contains information about the chr, the number of pairs. The second row and onwards list all the pairs
    loaded_segments_file = network_drawer.check_file_exist(
        separator="\t", header_value=None, skip_rows=1)

    # This just renames the first row to Pairs which can be used to identify the column later
    loaded_segments_file = loaded_segments_file.rename(columns={0: "Pairs"})

    # This just drops any empty rows in the dataframe
    if loaded_segments_file["Pairs"].isnull().any():
        loaded_segments_file = network_drawer.drop_empty_rows(
            loaded_segments_file)

    network_drawer_df = network_drawer.isolate_ids(
        loaded_segments_file, iid_list)

    carrier_in_network_dict = network_drawer.carriers_in_network(
        iid_list, network_drawer_df, ind_in_network_dict, variant_of_interest)

    network_drawer.draw_networks(network_drawer_df)

    return carrier_in_network_dict
