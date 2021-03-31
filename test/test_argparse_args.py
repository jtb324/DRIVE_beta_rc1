import argparse
import sys
sys.path.append("../drive")

import parser_arguments

def run():
    pass

def test_create_args():
    parser = parser_arguments.create_args(run)

    assert type(parser) == argparse.ArgumentParser, "the parser is not the correct type"