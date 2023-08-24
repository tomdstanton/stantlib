#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = 'Tom Stanton'
__title__ = 'Clean pairs'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Remove non-paired reads from paired fastq files'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import sys
import argparse
import gzip
from pathlib import Path
from typing import Generator, Iterable
import itertools
import datetime
import re

FASTQ_REGEX = re.compile('|'.join([
    r'_R[12]\.(f(?:ast)?q(?:\.gz)?)$',
    r'_R[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    r'_R[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    r'_[12]\.(f(?:ast)?q(?:\.gz)?)$',
    r'_[12]_[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$',
    r'_[12].[0-9]+?\.(f(?:ast)?q(?:\.gz)?)$']), re.IGNORECASE)

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
RED = '\033[31m'


class FastqFile:
    def __init__(self, path: Path | None, name: str | None = None, extension: str | None = None,
                 compressed: bool | None = None):
        self.path = path or Path()
        self.name = name or ''
        self.extension = extension or ''
        self.compressed = compressed or False

    @classmethod
    def from_path(cls, path: str):
        path = Path(path).absolute()
        if not path.exists():
            quit_with_error(f"File {path} does not exist")
        elif not path.is_file():
            quit_with_error(f"{path} is not a file")
        elif path.stat().st_size < 10:
            quit_with_error(f"{path} is less than 10 bytes")
        elif not (extension_match := FASTQ_REGEX.search(path.name)):
            quit_with_error(f"{path} is not a fastq file")
        else:
            return cls(path=path, name=path.name.replace(extension_match[0], ''), extension=extension_match[0],
                       compressed=path.suffix == '.gz')

    def __repr__(self):
        return self.name

    def __iter__(self) -> Generator:
        return parse_fastq(self.path, self.compressed)

    def __next__(self):
        return next(self.__iter__())


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


def parse_fastq(path: Path, compressed: bool = False) -> Generator[tuple[bytes, bytes, bytes, bytes], None, None]:
    """
    Open a fastq file and read 4 lines at a time using itertools zip_longest
    """
    open_func = gzip.open if compressed else open
    try:
        with open_func(path, 'rb') as stream:
            yield from itertools.zip_longest(*[iter(stream)] * 4, fillvalue=None)
    except Exception as e:
        quit_with_error(f"Error reading file {path}: {e}")


def get_headers(path: Path, compressed: bool = False) -> set[bytes]:
    """
    Open a fastq file and read the header line of each record using itertools islice
    """
    open_func = gzip.open if compressed else open
    try:
        with open_func(path, 'rb') as stream:
            return set(i.split(b' ', 1)[0] for i in itertools.islice(stream, 0, None, 4))
    except Exception as e:
        quit_with_error(f"Error reading file {path}: {e}")


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=__description__, formatter_class=argparse.RawTextHelpFormatter,
        add_help=False, usage='%(prog)s <fastq1> <fastq2> [options]',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')
    positionals = parser.add_argument_group(bold('Input'))
    positionals.add_argument('fastq1', help='Fastq file', type=lambda x: FastqFile.from_path(x))
    positionals.add_argument('fastq2', help='Fastq file', type=lambda x: FastqFile.from_path(x))
    options = parser.add_argument_group(bold('Options'))
    options.add_argument('-o', '--output', help='Output path for cleaned reads (default: %(default)s)',
                         metavar="", type=Path, default=Path.cwd())
    options.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    options.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    if len(a) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


def main():
    args = parse_args(sys.argv[1:])

    # Init the headers set with the headers from the first file
    log(f"Reading headers from {args.fastq1.path.name}")
    fastq1_headers = get_headers(args.fastq1.path, args.fastq1.compressed)
    headers_intersection = set()
    # We need to parse over the first file twice, but only parse over the second file once
    # We can use the headers from the second file to filter the set and filter the second file at the same time

    log(f"Filtering {args.fastq2.path.name}")
    with open(args.output / f"{args.fastq2}_cleaned{args.fastq2.extension}", 'wb') as output:
        for record in args.fastq2:
            if (header := record[0].split(b' ', 1)[0]) in fastq1_headers:
                output.write(b''.join(record))
                headers_intersection.add(header)

    if not headers_intersection:
        quit_with_error("No headers in common between the two files")

    log(f"Found {len(headers_intersection)} headers in common")

    log(f"Filtering {args.fastq1.path.name}")
    with open(args.output / f"{args.fastq1}_cleaned{args.fastq1.extension}", 'wb') as output:
        for record in args.fastq1:
            if (header := record[0].split(b' ', 1)[0]) in headers_intersection:
                output.write(b''.join(record))
                headers_intersection.remove(header)

    log("Done!")


if __name__ == '__main__':
    main()
