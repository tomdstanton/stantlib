#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = 'Tom Stanton'
__title__ = 'Genbank Slicer'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Filter a Genbank file by slicing or searching records and features'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import argparse
import sys
import re
import datetime

from Bio import SeqIO

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
RED = '\033[31m'


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


def check_slice(slice: str) -> tuple[int, int]:
    if not slice or not isinstance(slice, str) or len(slice) < 3:
        quit_with_error(f"Bad slice: {slice}")
    elif ":" not in slice:
        quit_with_error(f"Bad slice: {slice}")
    else:
        start, end = slice.split(":", 1)
        if not start.isdigit() or not end.isdigit():
            quit_with_error(f"Bad slice: {slice}")
        elif int(start) == int(end) or int(start) > int(end):
            quit_with_error(f"Bad slice: {slice}")
        else:
            return int(start), int(end)


def merge_ranges(ranges: list[tuple[int, int]], tolerance=0) -> list[tuple[int, int]]:
    """
    Merge overlapping ranges, ASSUMES EACH RANGE IS SORTED
    :param ranges: List of tuples of start and end positions
    :param tolerance: The number of bases to allow between alignments to be considered continuous
    :return: List of merged ranges
    """
    if ranges:
        sorted_ranges, merged_ranges = sorted(ranges, key=lambda x: x[0]), []
        # current_range = sorted(sorted_ranges[0])  # Sort the first range to ensure the start is before the end
        current_range = sorted_ranges[0]
        for start, end in sorted_ranges[1:]:
            # start, end = sorted((start, end))  # Sort the range to ensure the start is before the end
            if start - tolerance <= current_range[1]:
                # Ranges overlap, merge them
                current_range = (current_range[0], max(current_range[1], end))
            else:
                # No overlap, add the current range to the merged list and start a new range
                merged_ranges.append(current_range)
                current_range = (start, end)
        # Add the last range to the merged list
        merged_ranges.append(current_range)
        return merged_ranges
    else:
        return []


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=bold(__description__), add_help=False,
        usage='%(prog)s <genbank> [filter] [options] > out.{gbk,fa}',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')

    opts = parser.add_argument_group(bold('Input'))
    opts.add_argument('genbank', help='Genbank file or - for stdin', type=argparse.FileType('rt'))

    opts = parser.add_mutually_exclusive_group(required=True)
    opts.add_argument('-s', '--slice', metavar="", nargs="+", type=check_slice, default=[],
                      help='Slices each record by ranges given in format start:end')
    opts.add_argument('-f', '--feature', metavar="", type=lambda x: re.compile(x),
                      help="Filter by features matching regex pattern")
    opts.add_argument('-r', '--record', metavar="", type=lambda x: re.compile(x),
                      help="Filter by records matching regex pattern")

    opts = parser.add_argument_group(bold('Other options'))
    opts.add_argument('-o', '--outfmt', choices=['genbank', 'fasta'], default='genbank', metavar="",
                      help='Output either genbank or fasta format')
    opts.add_argument('-m', '--merge', type=int, default=0, metavar="",
                      help='Merge slice ranges within this distance, (default: 0)')
    opts.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    opts.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                      help='Show program version and exit')

    if len(a) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


if __name__ == '__main__':
    args = parse_args(sys.argv[1:])

    for record in SeqIO.parse(args.genbank, 'genbank'):
        feature_slices = [] + args.slice

        if args.feature:
            for feature in record.features:
                if args.feature.search(str(feature.qualifiers)):
                    feature_slices.append((feature.location.start, feature.location.end))

        if args.record:
            if args.record.search(str(record)):
                SeqIO.write(record, sys.stdout, args.outfmt)

        if feature_slices:
            for start, end in merge_ranges(feature_slices, args.merge):
                SeqIO.write(record[start:end], sys.stdout, args.outfmt)
