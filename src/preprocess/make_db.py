"""Make a glycan database byonic file from a csv file.

The csv file should have at least one column: "composition".
The "composition" column has glycan compositions in the condensed forms,
e.g. H5N4S2.

Usage:
    $ python make_db.py glycans.csv db.byonic
    #                   ~~~~~~~~~~~ ~~~~~~~~~
    #                      input     output
"""
import re
import csv


def convert(comp: str) -> str:
    """Convert composition from condensed format to byonic format.
    
    >>> convert('H3N4F2S1')
    'Hex(3)HexNAc(4)dHex(2)NeuAc(1)'
    >>> convert('H6N2')
    'Hex(6)HexNAc(2)'
    """
    comp = re.sub(r"H(\d+)", r"Hex(\1)", comp)
    comp = re.sub(r"N(\d+)", r"HexNAc(\1)", comp)
    comp = re.sub(r"F(\d+)", r"dHex(\1)", comp)
    comp = re.sub(r"S(\d+)", r"NeuAc(\1)", comp)
    return comp


input_fp = open(snakemake.input[0], "r", encoding="utf8")
output_fp = open(snakemake.output[0], "w", encoding="utf8")
with input_fp, output_fp:
    reader = csv.reader(input_fp)
    next(reader)
    for comp, _ in reader:
        byonic_comp = convert(comp)
        row = byonic_comp + " %\n"
        output_fp.write(row)
