#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = 'Tom Stanton'
__title__ = 'GC Filter'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Filter fastq file by GC content'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import sys
import argparse
import gzip
from typing import Generator, Iterable
import itertools
from pathlib import Path
import datetime
from tempfile import NamedTemporaryFile


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


def decompress_byte_stream(byte_stream: bytes) -> bytes:
    if not byte_stream:
        quit_with_error("Nothing to decompress")
    try:
        decompressed_data = gzip.decompress(byte_stream)
        log("Decompressed input")
        return decompressed_data
    except gzip.BadGzipFile:
        log("Not a gzip file, assuming uncompressed")
        return byte_stream
    except OSError:
        quit_with_error("Error decompressing file")


def read_file(file: str, compressed: bool | str = "unknown") -> bytes:
    if file == '-':
        log("Reading from stdin")
        # sys.stdin.seek(0)  # We are potentially reading from buffer twice so make sure it's at the beginning
        if compressed == "unknown" or True:
            return decompress_byte_stream(sys.stdin.buffer.read())
        else:
            return sys.stdin.buffer.read()
    else:
        if (file := Path(file)).is_file() and file.stat().st_size > 0:
            log(f"Reading from {file}")
            compressed = file.suffix == ".gz" if compressed == "unknown" else compressed
            if compressed:
                return decompress_byte_stream(file.read_bytes())
            else:
                return file.read_bytes()
        else:
            quit_with_error(f"Error reading file {file}")


def fastq_iter(byte_stream: bytes) -> Generator:
    """Iterate over fastq file using itertools to yield a split on every 4th line"""
    if byte_stream.startswith(b'@'):
        return itertools.zip_longest(*[iter(byte_stream.splitlines())] * 4, fillvalue=None)
    else:
        quit_with_error("File does not appear to be fastq")


def parse_fastq(file: str, seq_only: bool = False, **kwargs) -> Generator:
    return (r[1] if seq_only else r for r in fastq_iter(read_file(file, **kwargs)) if r)


def count_gc(seq: bytes) -> float:
    return (seq.count(b'G') + seq.count(b'C')) / len(seq)


def get_stats(floats: Iterable[float]) -> tuple[int, float, float]:
    """Returns mean and standard deviation of an iterable of floats
    Tries to iterate over the floats once"""
    log("Calculating GC content stats")
    count, mean, M2 = 0, 0, 0
    for x in floats:
        count += 1
        delta = x - mean
        mean += delta / count
        M2 += delta * (x - mean)
    if count < 2:
        quit_with_error("Not enough reads to calculate stats")
    else:
        return count, mean, (M2 / (count - 1)) ** 0.5


def draw_hist(data: Iterable[float], num_bins: int = 100, value_range: tuple[float, float] = (0, 1)):
    """
    Draws a simple histogram in the terminal
    """
    # Create histogram bins
    bin_width = (value_range[1] - value_range[0]) / num_bins
    bins = [value_range[0] + i * bin_width for i in range(num_bins + 1)]
    counts = [0] * num_bins
    for value in data:
        bin_index = min(int((value - value_range[0]) / bin_width), num_bins - 1)
        counts[bin_index] += 1
    max_count = max(counts)
    plot = []
    for c in counts:
        if max_count == 0:
            normalized_count = 0
        else:
            normalized_count = int(c / max_count * 40)  # Scale to 40 characters width
        plot.append("*" * normalized_count)
    print("GC Bins:")
    for i in range(num_bins):
        if counts[i] > 0:
            bin_label = f"{bins[i]:.2f}-{bins[i + 1]:.2f}"
            print(f"{bin_label}: {plot[i]}")


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=__description__, formatter_class=argparse.RawTextHelpFormatter,
        add_help=False, usage='%(prog)s <fastq> [options]',
        epilog='If --lower and --upper are not specified, stats will be used to determine the cutoffs')
    input = parser.add_argument_group(bold('Input'))
    input.add_argument('fastq', help='Path to fastq(.gz) file or - for stdin', metavar="<fastq>",
                       # nargs='?', default='-'
                       )
    options = parser.add_argument_group(bold('Options'))
    options.add_argument('-u', '--upper', type=float, default=1, help='Upper GC content cutoff as float',
                         metavar="1.0")
    options.add_argument('-l', '--lower', type=float, default=0, help='Lower GC content cutoff as float',
                         metavar="0.0")
    options.add_argument('-s', '--stats', action='store_true', help='Print stats and exit')
    options.add_argument('-r', '--reverse', action='store_true', help='Reverse filter')
    options.add_argument('-p', '--print', action='store_true', help='Print the GC of each read and exit')
    options.add_argument('-d', '--draw', nargs="?", const=100, default=None, type=int, metavar='100',
                         help='Draw a histogram of read GC and exit\nProviding an int will change binwidth')
    options.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    options.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    if len(a) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


# def main():
if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    # args = parse_args(['/Users/tom/Bioinformatics/stantlib/stantlib/test_data/ovc1_1.fastq'])

    if args.print:
        for _, seq, _, _, in parse_fastq(args.fastq):
            print(count_gc(seq))
        log("Done!")
        sys.exit(0)

    if args.draw:
        draw_hist((count_gc(seq) for _, seq, _, _, in parse_fastq(args.fastq)), num_bins=args.draw)
        log("Done!")
        sys.exit(0)

    if args.stats:
        count, mean, std = get_stats(count_gc(seq) for _, seq, _, _, in parse_fastq(args.fastq))
        print(f"Reads: {count}\nMean GC: {mean * 100:.2f}%\nStd: {std * 100:.2f}%")
        log("Done!")
        sys.exit(0)

    if args.lower == 0 and args.upper == 1:
        tempfile = NamedTemporaryFile()
        gc = []
        with open(tempfile.name, "wb") as f:
            for record in parse_fastq(args.fastq):
                f.write(b'\n'.join(record) + b'\n')
                gc.append(count_gc(record[1]))
        stats = get_stats(gc)
        args.lower = stats[1] - stats[2]
        args.upper = stats[1] + stats[2]
        args.fastq = tempfile.name

    if args.lower > args.upper:
        quit_with_error("Lower cutoff is greater than upper cutoff")

    if not args.reverse:
        log(f"Filtering reads with GC content between {args.lower * 100:.2f}% and {args.upper * 100:.2f}%")
    else:
        log(f"Filtering reads with GC content < {args.lower * 100:.2f}% and > {args.upper * 100:.2f}%")

    for record in parse_fastq(args.fastq):
        if args.lower <= count_gc(record[1]) <= args.upper and not args.reverse:
            sys.stdout.buffer.write(b'\n'.join(record) + b'\n')
        elif not args.lower <= count_gc(record[1]) <= args.upper and args.reverse:
            sys.stdout.buffer.write(b'\n'.join(record) + b'\n')

    log("Done!")
    sys.exit(0)
