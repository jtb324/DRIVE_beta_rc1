import logging

from network_creator_class import Network_Img_Maker


def create_networks(segments_file, output_path):

    logger = logging.getLogger(output_path+'/single_variant_analysis.log')

    network_drawer = Network_Img_Maker(segments_file, output_path, logger)

    loaded_segments_file = network_drawer.check_file_exist(
        separator=",", header_value=None, skip_rows=5)

    network_drawer_df = network_drawer.isolate_ids(loaded_segments_file)

    network_drawer.draw_networks(network_drawer_df)
