#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = 'Tom Stanton'
__title__ = 'genbank2seq'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Extract sequences from Genbank records'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import sys
import datetime
import argparse
from Bio import SeqIO

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
RED = '\033[31m'
YELLOW = '\033[93m'


# Functions --------------------------------------------------------------------
def bold(text: str):
    return BOLD + text + END_FORMATTING


def bold_red(text: str):
    return RED + BOLD + text + END_FORMATTING


def bold_yellow(text: str):
    return YELLOW + BOLD + text + END_FORMATTING


def warning(text: str):
    log(bold_yellow(f"WARNING: {text}"))


def error(text: str):
    log(bold_red(f"ERROR: {text}"))


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def log(message: str = '', end: str = '\n', sep: str = ' ', flush: bool = True, file=sys.stderr):
    print(f"{get_timestamp()}: {message}", file=file, flush=flush, end=end, sep=sep)


def quit_with_error(message: str):
    error(message)
    sys.exit(1)


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=bold(__description__), add_help=False, usage=f'%(prog)s <genbank> [options] > output.fa',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')

    positionals = parser.add_argument_group(bold('Input'))
    positionals.add_argument('genbank', help='Genbank file or - for stdin', type=argparse.FileType('rt'))

    record_options = parser.add_argument_group(bold("Records"))
    record_options.add_argument('-r', '--record', default='locus', choices=['locus', 'feature'],
                                help='Extract the sequence from each locus or feature (default: locus)', metavar='')
    record_options.add_argument('-f', '--feature_type', default='CDS', metavar='',
                                help='If --record=feature, which feature type to extract (default: CDS)')
    record_options.add_argument('-i', '--identifier', default='locus_tag', metavar='',
                                help='If --record=feature, which qualifier to use as the >id (default: locus_tag)')

    translation_options = parser.add_argument_group(bold('Protein'))
    translation_options.add_argument('-p', '--protein', action='store_true',
                                     help='If --record=feature, extract the feature translation')
    translation_options.add_argument('-t', '--translate', action='store_true',
                                     help='Translate on-the-fly, turns on --protein')
    translation_options.add_argument('-S', '--to_stop', action='store_true',
                                     help='If --translate, translate to the first stop codon')
    translation_options.add_argument('-s', '--stop_symbol', default='*', metavar='',
                                     help='If --translate, use this symbol for stop codons (default: *)')
    translation_options.add_argument('-T', '--table', default=1, type=int, metavar='',
                                        help='If --translate, use this translation table (default: 1)')

    other_options = parser.add_argument_group(bold('Other'))
    other_options.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other_options.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                               help='Show version number and exit')

    if len(a) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    if args.translate:
        args.protein = True
    if args.record == 'locus' and args.protein:
        quit_with_error('Cannot extract protein sequence from locus')

    for record in SeqIO.parse(args.genbank, 'genbank'):
        if args.record == 'locus':
            SeqIO.write(record, sys.stdout, 'fasta')
        else:
            for feature in record.features:
                if feature.type == args.feature_type:
                    if args.identifier not in feature.qualifiers:
                        warning(f'{args.identifier} qualifier not found for {feature}')
                        continue
                    else:
                        identifier = feature.qualifiers[args.identifier][0]

                    if args.protein:
                        if args.translate:
                            seq = feature.extract(record.seq).translate(
                                to_stop=args.to_stop, table=args.table, stop_symbol=args.stop_symbol)
                        else:
                            if 'translation' in feature.qualifiers:
                                seq = feature.qualifiers['translation'][0]
                            else:
                                warning(f'No translation found for {identifier}')
                                continue
                    else:
                        seq = feature.extract(record.seq)

                    sys.stdout.write(f">{identifier}\n{seq}\n")
