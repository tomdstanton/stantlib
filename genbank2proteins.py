#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
genbank2proteins
-----------------------------------------------------------------------------------------
Program to extract protein translations from a Genbank file.
-----------------------------------------------------------------------------------------
Requires:
- Python >=3.9
- biopython
-----------------------------------------------------------------------------------------
Tom Stanton, 2023
"""

from __future__ import annotations

import sys, argparse

from Bio import SeqIO

DESCRIPTION = """
Extract translations from a Genbank file
-----------------------------------------------------------------------------------------
    Takes a Genbank records, extracts all protein translations from each
    feature using the translation attribute or on the fly, and returns
    amino acid sequences in fasta format.
-----------------------------------------------------------------------------------------
"""


def parse_args(a):
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter,
                                     usage="python3 genbank2protein.py <in.gbk> <out.faa>")
    parser.add_argument("-i", "--input", default=sys.stdin, type=argparse.FileType('r'),
                        help="Input file [default: stdin]\n"
                             "The translation in the Genbank features will be used\n"
                             "unless --translate is used, in which case the\n"
                             "DNA sequence will be translated.")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType('w', encoding='UTF-8'),
                        help="Output path [default: stdout]\n"
                             "If a path is given, each translation will be written to a separate\n"
                             "file in the directory named after the locus tag.")
    parser.add_argument("--translate", action="store_true", default=False,
                        help="DNA features will be extracted and translated with BioPython.\n"
                             "If False, the translation in the Genbank features will be used.")
    parser.add_argument("--table", default=11, help="Translation table for BioPython translate", metavar='[11]')
    parser.add_argument("--to_stop", action="store_true", default=False,
                        help="Stop translation at first stop codon for BioPython translate")
    parser.add_argument("--stop_symbol", type=str, default="*", metavar="[*]",
                        help="Stop symbol for translation")
    parser.add_argument("--feature_type", type=str, default="CDS", metavar="[CDS]",
                        help="Genbank feature type to extract")
    return parser.parse_args(a)


def log(message: str, flush: bool = True, out=sys.stderr, end='\n'):
    print(message, flush=flush, file=out, end=end)


def warning(message: str, **kwargs):
    message = f"\033[1;33m{message}\033[0m"  # Wrap message in bold yellow text
    log(f"WARNING: {message}", **kwargs)


def error(message: str, **kwargs):
    message = f"\033[1;31m{message}\033[0m"  # Wrap message in bold red text
    log(f"ERROR: {message}", **kwargs)


def quit_with_error(message: str, **kwargs):
    error(message, **kwargs)
    sys.exit(1)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    for record in SeqIO.parse(args.input, 'genbank'):
        for feature in record.features:
            if feature.type == args.feature_type:
                print(feature.qualifiers)
                identifier = feature.qualifiers['locus_tag'][0]
                if args.translate:
                    seq = feature.extract(record.seq).translate(
                        to_stop=args.to_stop, table=args.table, stop_symbol=args.stop_symbol)
                elif 'translation' in feature.qualifiers:
                    seq = feature.qualifiers['translation'][0]
                else:
                    warning(f'No translation found for {identifier}')
                    continue
                args.output.write(f'>{identifier}\n{seq}\n')
