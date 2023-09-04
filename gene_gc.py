#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Tom Stanton'
__title__ = 'Gene GC'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Calculate GC content of genes in a genbank file'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import sys
import argparse
import datetime

from Bio import SeqIO

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
RED = '\033[31m'


# Functions -------------------------------------------------------------------
def bold(text: str):
    return BOLD + text + END_FORMATTING


def bold_red(text: str):
    return RED + BOLD + text + END_FORMATTING


def error(text: str):
    log(bold_red(f"ERROR: {text}"))


def get_timestamp():
    return '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())


def log(message: str = '', end: str = '\n', sep: str = ' ', flush: bool = True, file=sys.stderr):
    print(f"{get_timestamp()}: {message}", file=file, flush=flush, end=end, sep=sep)


def quit_with_error(message: str):
    error(message)
    sys.exit(1)


def get_gc(seq: str) -> tuple[int, int]:
    seq = seq.upper()
    return seq.count('G'), seq.count('C')


def get_gc_sum(seq: str) -> int:
    return sum(get_gc(seq))


def get_gc_frac(seq: str) -> float:
    if len(seq) == 0:
        return 0
    return get_gc_sum(seq) / len(seq)


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=bold(__description__), formatter_class=argparse.RawTextHelpFormatter,
        add_help=False, usage='%(prog)s <genbank> [options] > out.tab',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')
    opts = parser.add_argument_group(bold('Input'))
    opts.add_argument('genbank',
                      help='Path to genbank file or ', type=argparse.FileType('rt'))
    opts = parser.add_argument_group(bold('Options'))
    opts.add_argument('-f', '--feature_type', metavar="", default="CDS", type=str,
                      help='Feature type to report (default: %(default)s)')
    opts.add_argument('-s', '--suppress_header', action='store_true',
                      help='Suppress header in tab-delimited output')
    opts.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    opts.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                      help='Show program version and exit')

    if len(a) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    header_printed = args.suppress_header
    for r in SeqIO.parse(args.genbank, 'genbank'):
        bg_gc = get_gc_frac(r.seq)

        if not header_printed:
            sys.stdout.write(
                "reference\treference_gc\tfeature_type\tfeature_start\tfeature_end\t"
                "feature_name\tfeature_locus_tag\tfeature_product\tgc\tgc_diff\n"
            )
        for f in r.features:
            if f.type == args.feature_type:
                gc = get_gc_frac(f.extract(r).seq)
                name = f.qualifiers.get('gene', [''])[0]
                locus_tag = f.qualifiers.get('locus_tag', [''])[0]
                desc = f.qualifiers.get('product', [''])[0]
                sys.stdout.write(
                    f"{r.id}\t{bg_gc:.3f}\t{args.feature_type}\t{name}\t{locus_tag}\t{desc}\t{f.location.start}\t"
                    f"{f.location.end}\t{gc:.3f}\t{gc - bg_gc:.3f}\n"
                )

    log("Done!")
