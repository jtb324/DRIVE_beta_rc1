#!/usr/bin/env python
import typer

from carrier_identification import carrier_id_main
from segment_analysis import segment_analysis_main

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
