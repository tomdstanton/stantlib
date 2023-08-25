#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = 'Tom Stanton'
__title__ = 'xccessions'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Convert a list of fixed-length accessions to a regex'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import re
import sys
import datetime
from typing import Iterable
import argparse

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


def check_equal_length(strings: Iterable[str]) -> Iterable[str]:
    if len(set(map(len, strings))) == 1:
        return strings
    else:
        quit_with_error("Accessions are not all the same length")


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=__description__, formatter_class=argparse.RawTextHelpFormatter,
        add_help=False, usage='%(prog)s <accessions> [options]',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')

    positionals = parser.add_argument_group(bold('Input'), "Breaks are relative to the reference")
    positionals.add_argument('accessions', nargs='+', help='List of accessions to convert to a regex',
                             type=check_equal_length)
    other_options = parser.add_argument_group(bold("Other options"))
    other_options.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other_options.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                               help='Show version number and exit')
    if len(a) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


def main():
    args = parse_args(sys.argv[1:])
    regex = ''
    for i in range(len(args.accessions[0])):
        chars = set(map(lambda x: x[i], args.accessions))
        if len(chars) == 1:
            regex += chars.pop()
        else:
            regex += f"[{''.join(chars)}]"
    try:
        r = re.compile(regex)
        matches = list(filter(r.match, args.accessions))  # Use the pattern to match all the input accessions
        if len(matches) == len(args.accessions):
            print(regex)
            log("Done!")
        else:
            quit_with_error(f"Regex pattern does not match all accessions: {regex}")
    except re.error:
        quit_with_error("Failed to compile regex: {regex_pattern}")


if __name__ == '__main__':
    main()
