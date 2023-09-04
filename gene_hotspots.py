#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Tom Stanton'
__title__ = 'Gene Hotspots'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Find DNA metric hotspots and report overlapping features'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import sys
import argparse
from typing import Generator, Iterable
import datetime
import collections
import math

from Bio import SeqIO

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
RED = '\033[31m'

# See: https://doc.ugene.net/wiki/display/UM38/DNA+Flexibility
DINUCL_ANGLES = {"AA": 7.6, "CA": 14.6, "AC": 10.9, "CC": 7.2, "AG": 8.8, "CG": 11.1, "AT": 12.5, "CT": 8.8,
                 "GA": 8.2, "TA": 25, "GC": 8.9, "TC": 8.2, "GG": 7.2, "TG": 14.6, "GT": 10.9, "TT": 7.6}


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


def get_windows(seq: list | str | bytes, window_size: int, step_size: int
                ) -> Generator[tuple[int, int, list | str | bytes], None, None]:
    """
    Returns a generator of windows of a given size and step size
    """
    for start in range(0, len(seq), step_size):
        yield start, (stop := start + window_size), seq[start:stop]


def get_gc(seq: str) -> tuple[int, int]:
    seq = seq.upper()
    return seq.count('G'), seq.count('C')


def get_gc_sum(seq: str) -> int:
    return sum(get_gc(seq))


def get_gc_frac(seq: str) -> float:
    if len(seq) == 0:
        return 0
    return get_gc_sum(seq) / len(seq)


def get_gc_skew(seq: str) -> float:
    g, c = get_gc(seq)
    if g + c == 0:
        return 0
    return (g - c) / (g + c)


def get_flexibility(seq: str) -> float:
    """
    Calculate the DNA flexibility of a sequence using:
    (average seq threshold) = (sum of flexibility angles in the seq) / (seq len - 1)
    """

    return sum(DINUCL_ANGLES[s] for s in get_windows(seq.upper(), 2, 1)) / (len(seq) - 1)


def get_stats(floats: Iterable[float | int]) -> tuple[float, float]:
    """Returns mean and standard deviation of an iterable of floats or ints in a single pass."""
    floats = list(floats)
    if not floats:
        quit_with_error("No data to calculate stats")
    mean = sum(floats) / len(floats)
    std = math.sqrt(sum((x - mean) ** 2 for x in floats) / len(floats))
    return mean, std


def get_entropy(s):
    return sum(-p_x * math.log(p_x, 2) for p_x in [n_x / len(s) for x, n_x in collections.Counter(s).items()])


def process_subseq(subseq: str, analysis: str) -> float:
    if analysis == "gc_frac":
        return get_gc_frac(subseq)
    elif analysis == "gc_sum":
        return get_gc_sum(subseq)
    elif analysis == "gc_skew":
        return get_gc_skew(subseq)
    elif analysis == "entropy":
        return get_entropy(subseq)
    elif analysis == "flexibility":
        return get_flexibility(subseq)
    else:
        quit_with_error(f"Unknown analysis: {analysis}")


def find_hotspots(windows: dict[tuple[int, int], float], zscore: float) -> Generator[tuple[tuple[int, int], float], None, None]:
    """
    Given a dictionary of windows where keys represent the start and end of the window and the value is a metric
    such as gc content or entropy, find windows where the metric is greater than the zscore threshold.
    """
    # First we need to calculate the mean and standard deviation of the metric
    mean, std = get_stats(windows.values())
    upper_threshold = mean + (std * zscore)
    lower_threshold = mean - (std * zscore)
    log(f"mean: {mean:.2f}\tstd: {std:.2f}\tupper_threshold: {upper_threshold:.2f}\tlower_threshold: {lower_threshold:.2f}")
    for k, v in windows.items():
        if v > upper_threshold or v < lower_threshold:
            yield k, v


def merge_ranges(ranges: Iterable[tuple[int, int]], tolerance=0) -> list[tuple[int, int]]:
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


def range_overlap(range1: tuple[int, int], range2: tuple[int, int], skip_sort: bool = False) -> int:
    """
    Returns the overlap between two ranges
    :param range1: Tuple of start and end positions
    :param range2: Tuple of start and end positions
    :param skip_sort: Skip sorting each range before calculating the overlap
    :return: Integer of overlap
    """
    start1, end1 = range1 if skip_sort else sorted(range1)
    start2, end2 = range2 if skip_sort else sorted(range2)
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=bold(__description__), formatter_class=argparse.RawTextHelpFormatter,
        add_help=False, usage='%(prog)s <genbank> -a analysis [options] > out.tab',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')

    opts = parser.add_argument_group(bold('Input'))
    opts.add_argument('genbank', help='Genbank file or - for stdin', type=argparse.FileType('rt'))

    opts = parser.add_argument_group(bold('Options'))
    opts.add_argument('-a', '--analysis', metavar="", default="gc_frac", type=str,
                      choices=['gc_frac', 'gc_sum', 'gc_skew', 'entropy', 'flexibility'],
                      help="Analysis to perform (default: %(default)s)\n - gc_frac: GC fraction\n - gc_sum: GC sum\n"
                           " - gc_skew: GC skew\n - entropy: Shannon entropy\n - flexibility: DNA flexibility")
    opts.add_argument('-w', '--window-size', metavar="", default=500, type=int,
                      help='Window size (default: %(default)s)')
    opts.add_argument('-s', '--step-size', metavar="", default=250, type=int,
                      help='Step size (default: %(default)s)')
    opts.add_argument('-z', '--zscore', metavar="", default=3, type=float,
                      help='Z-score threshold (default: %(default)s)')
    opts.add_argument('-m', '--merge', metavar="", default=100, type=int,
                      help='Range merge tolerance (default: %(default)s)')
    opts.add_argument('-f', '--feature_type', metavar="", default="CDS", type=str,
                      help='Feature type to report (default: %(default)s)')
    opts.add_argument('--suppress_header', action='store_true',
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
        # Get window data
        data = {}
        for start, stop, subseq in get_windows(r.seq, args.window_size, args.step_size):
            data[(start, stop)] = process_subseq(subseq, args.analysis)
        # Find hotspots
        log(f"Finding {args.analysis} hotspots in {r.id} using a z-score threshold of {args.zscore}")
        hotspots = dict(find_hotspots(data, zscore=args.zscore))
        if hotspots:
            if not header_printed:
                sys.stdout.write(
                   "reference\tstart\tstop\tanalysis\tvalue\tfeature_type\tfeature_start\tfeature_end\tfeature_name\t"
                   "feature_locus_tag\tfeature_description\n"
                )
                header_printed = True

            merged_hotspots = {}  # Merge overlapping hotspots
            for start, stop in merge_ranges(hotspots.keys(), args.merge):
                merged_hotspots[(start, stop)] = process_subseq(r.seq[start:stop], args.analysis)

            log(f"Found {len(merged_hotspots)} hotspots")

            # See if and features overlap with the GC hotspots
            for h_range, v in merged_hotspots.items():
                for f in r.features:
                    if f.type == args.feature_type and range_overlap(h_range, (f.location.start, f.location.end)) > 0:
                        name = f.qualifiers.get('gene', [''])[0]
                        locus_tag = f.qualifiers.get('locus_tag', [''])[0]
                        desc = f.qualifiers.get('product', [''])[0]
                        sys.stdout.write(
                            f"{r.id}\t{h_range[0]}\t{h_range[1]}\t{args.analysis}_hotspot\t{v:.2f}\t{f.type}\t"
                            f"{f.location.start}\t{f.location.end}\t{name}\t{locus_tag}\t{desc}\n"
                        )
    log("Done!")
