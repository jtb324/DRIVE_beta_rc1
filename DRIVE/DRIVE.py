#!/usr/bin/env python
import typer
from dotenv import load_dotenv

from carrier_identification import carrier_id_main
from segment_analysis import segment_analysis_main

# load in the environment variables from the drive.env file
load_dotenv("../drive.env")

# creating a cli object
cli = typer.Typer(
    add_completion=False,
    help="Tool identifies networks of individuals who share an IBD segment around a genomic region of interest",
)
# create a subcommand for the carrier identification step
cli.add_typer(carrier_id_main.carrier_id_app, name="carrier_identification")

cli.add_typer(segment_analysis_main.segment_analysis_app, name="segment_analysis")


if __name__ == "__main__":
    cli()
