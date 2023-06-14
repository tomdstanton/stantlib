#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Tom Stanton'
__title__ = 'FSF'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Flatten, simplify, and filter a fasta file'
__license__ = 'gpl-3.0'
__version__ = '0.0.1'

# ................ Python Imports ................ #
import argparse
import re
import sys
from argparse import ArgumentParser
from pathlib import Path
from re import compile

# ................ Globals ................ #
REGEX = compile('[^a-zA-Z]')


class Rec:
    def __init__(self, header, description, sequence, args):
        self.sequence = REGEX.sub("", sequence) if args.strip else sequence
        self.header = f"{header} {description}" if args.description else header
        if args.case == 'upper':
            self.sequence = self.sequence.upper()
        elif args.case == 'lower':
            self.sequence = self.sequence.lower()
    
    def print_out(self):
        print(f">{self.header}")
        print(self.sequence)



def seq_to_rec(seq: str, args: argparse.Namespace) -> dict:
    seqs = {}
    for record in seq.split('\n>'):
        if record:
            lines = record.strip().split('\n')
            if lines:
                header = lines[0].split(' ')[0].replace(">", "")
                if header not in args.remove:
                    seqs[header] = Rec(
                        header=header,
                        description=lines[0].split(' ', 1)[1] if ' ' in lines[0] else '',
                        sequence=''.join(lines[1:]),
                        args=args
                        )
    return seqs


#
# def seq_to_dict(seq: str) -> dict:
#     return {r[0]: ''.join(r[1:]) for record in seq.split('>')[1:] if (r := record.splitlines())}


# def flatten_seq(seq: str | bytes) -> str | bytes:
#     for record in seq.split(b'>')[1:] if isinstance(seq, bytes) else seq.split('>')[1:]:
#         if lines := record.splitlines():
#             yield b'>' + b'\n'.join([lines[0], b''.join(lines[1:])]) \
#                 if isinstance(seq, bytes) else '>' + '\n'.join([lines[0], ''.join(lines[1:])])


# def dict_to_seq(seq: dict) -> str:
#     return '\n'.join([f'>{k}\n{v}' for k, v in seq.items()])



def parse_args(args):
    parser = ArgumentParser(description=f"\n{__title__} {__version__}: {__description__}", 
    usage=f"python3 {__file__} <fasta> -k $(cat accessions.list) -c upper ")
    parser.add_argument('fasta', type=lambda i: check_input(parser, i), default='-', help="Fasta file / stream")
    parser.add_argument('-d', '--description', action='store_true',  help="Include description in header")
    parser.add_argument('-s', '--strip', action='store_true', help="Strip invalid characters from sequence")
    parser.add_argument('-c', '--case', choices=['upper', 'lower'], help="Force sequence to specific case")
    parser.add_argument('-k', '--keep', nargs='+', default=[], type=str, help="List of accessions to keep")
    parser.add_argument('-r', '--remove', nargs='+', default=[], type=str, help="List of accessions to remove")
    return parser.parse_args(args)


def check_input(parser, path: str):
    if path == '-':
        return sys.stdin
    else:
        path = Path(path).absolute()
        if not path.stat().st_size:
            parser.error(f"The file {path} does not exist!")
        else:
            return open(path, 'r')


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    seqs = seq_to_rec(args.fasta.read(), args)
    args.fasta.close()
    if not seqs:
        print("No sequences found in file!", file=sys.stderr)
        sys.exit(1)

    if not args.keep:
        args.keep = list(seqs.keys())

    for acc in args.keep:
        seqs[acc].print_out()

    sys.exit(0)

