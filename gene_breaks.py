#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = 'Tom Stanton'
__title__ = 'Gene Breaks'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Find where an assembly has broken and report overlapping features'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import sys
import os
import datetime
from pathlib import Path
from typing import Iterable, Generator
from operator import attrgetter, le, ge
from time import time
import subprocess
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


def check_programs(progs: list[str], verbose: bool = False):
    """Check if programs are installed and executable"""
    bins = {  # Adapted from: https://unix.stackexchange.com/a/261971/375975
        f: Path(f'{p}/{f}') for p in filter(
            os.path.isdir, os.environ["PATH"].split(os.path.pathsep)
        ) for f in os.listdir(p) if os.access(f'{p}/{f}', os.X_OK)
    }
    for program in progs:
        if program in bins.keys():
            if verbose:
                log(f'{program}: in path {bins[program]}')
        else:
            quit_with_error(f'{program} not found')


def check_file(file: str) -> Path:
    if not (file := Path(file)).is_file():
        quit_with_error(f'{file.name} is not a file')
    elif file.stat().st_size == 0:
        quit_with_error(f'{file.name} is empty')
    else:
        return file.absolute()


def run_command(command: str, pipe: bytes = b'', shell: bool = False, cmd_split: str = ' ', verbose: bool = False,
                string_out: bool = False):
    """Run a command and return the output"""
    start = time()
    if verbose:
        log(f'Running Command\n\tPipe: {"True" if pipe else "False"}\n'
            f'\tShell: {"True" if shell else "False"}\n\tCommand: {command}')
        log(f'\r\tRunning...', flush=False, end=' ')
    if pipe:
        with subprocess.Popen(command.split(cmd_split), shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              stdin=subprocess.PIPE) as child:
            out, err = child.communicate(input=pipe)
    else:
        with subprocess.Popen(command.split(cmd_split), shell=shell, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE) as child:
            out, err = child.communicate()
    if verbose:
        log(f'Done! ({time() - start:.2f}s)')
        if err and not out:
            warning(err.decode())

    return out.decode() if string_out else out


def decode_line(line: str | bytes):
    return line.decode().strip() if isinstance(line, bytes) else line.strip()


def get_best_alignments(alignments: Iterable['Alignment'], group: str = "query_name", metric: str = "alignment_score",
                        sort_large_to_small: bool = True, threshold: float = 0.0) -> Generator['Alignment', None, None]:
    """
    Get the best alignments for each group and can additionally filter by calling filter_alignments internally
    """
    current_group = None
    for alignment in sorted(
            filter_alignments(alignments, metric, threshold, sort_large_to_small),
            key=lambda x: attrgetter(group, metric)(x), reverse=sort_large_to_small
    ):
        if current_group != getattr(alignment, group):
            current_group = getattr(alignment, group)
            yield alignment


def filter_alignments(alignments: Iterable['Alignment'], metric: str = "alignment_score", threshold: float = 0.0,
                      above_threshold: bool = True) -> Generator['Alignment', None, None]:
    """
    Filter alignments based on a metric and threshold
    """
    op = ge if above_threshold else le
    for alignment in alignments:
        assert hasattr(alignment, metric), f"{alignment} does not have attribute {metric}"
        assert isinstance(getattr(alignment, metric), (int, float)), f"{alignment} attribute {metric} is not a number"
        if magic(getattr(alignment, metric), threshold, op):
            yield alignment


def magic(left: int | float, right: int | float, op):
    return op(left, right)


def cull_conflicting_alignments(alignment_to_keep: 'Alignment', alignments: Iterable['Alignment']
                                ) -> Generator['Alignment', None, None]:
    """
    Returns a generator of alignments that do not conflict with the alignment_to_keep.
    param alignment_to_keep: Alignment object
    param alignments: Iterable of Alignment objects
    return: Generator of Alignment objects
    """
    for x in alignments:
        if not x.conflicts(alignment_to_keep):
            yield x


def cull_all_conflicting_alignments(alignments: Iterable['Alignment'], sort_by: str = 'alignment_score',
                                    sort_large_to_small: bool = True) -> list['Alignment']:
    kept_alignments = []
    sorted_alignments = sorted(list(alignments), key=lambda x: getattr(x, sort_by), reverse=sort_large_to_small)
    while sorted_alignments:
        kept_alignments.append(sorted_alignments.pop(0))
        sorted_alignments = list(cull_conflicting_alignments(kept_alignments[-1], sorted_alignments))
    return kept_alignments


def minimap2(
        query: Path | list[Path] | str | bytes, target: Path | str | bytes, preset: str = "", threads: int = 1,
        extra_args: str = "", verbose: bool = False) -> Generator['Alignment', None, None] | None:
    pipe = None
    query_stdin = isinstance(query, str) or isinstance(query, bytes)
    target_stdin = isinstance(target, str) or isinstance(target, bytes)
    if query_stdin and target_stdin:
        quit_with_error(f"Query and target cannot both be stdin")
    command = f"minimap2 -t {threads} "
    if preset:
        command += f"-x {preset} "
    if extra_args:
        command += f"{extra_args} "
    if query_stdin:
        command += f"{target} -"
        pipe = query.encode() if isinstance(query, str) else query
    elif query_stdin or target_stdin:
        command += f"- {' '.join(str(j) for j in query) if isinstance(query, list) else str(query)}"
        pipe = target.encode() if isinstance(target, str) else target
    else:
        command += f"{target} {' '.join(str(j) for j in query) if isinstance(query, list) else str(query)}"
    out = run_command(command, verbose=verbose, pipe=pipe) if pipe else run_command(command, verbose=verbose)
    return (Alignment.from_paf_line(i) for i in out.splitlines())


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


# Classes ----------------------------------------------------------------------
class AlignmentError(Exception):
    pass


class Alignment:
    def __init__(self, line: str | None = None, query_name: str | None = None, query_length: int | None = None,
                 query_start: int | None = None, query_end: int | None = None, strand: str | None = None,
                 target_name: str | None = None, target_length: int | None = None, target_start: int | None = None,
                 target_end: int | None = None, matching_bases: int | None = None, num_bases: int | None = None,
                 percent_identity: float | None = None, percent_query_coverage: float | None = None,
                 percent_target_coverage: float | None = None, alignment_score: float | None = None):

        self.query_name = query_name or ''
        self.query_length = query_length or 0
        self.query_start = query_start or 0
        self.query_end = query_end or 0
        self.strand = strand or ''
        self.target_name = target_name or ''
        self.target_length = target_length or 0
        self.target_start = target_start or 0
        self.target_end = target_end or 0
        self.matching_bases = matching_bases or 0
        self.num_bases = num_bases or 0
        self.percent_identity = percent_identity or 0
        self.percent_query_coverage = percent_query_coverage or 0
        self.percent_target_coverage = percent_target_coverage or 0
        self.alignment_score = alignment_score or 0

    @classmethod
    def from_paf_line(cls, line: str | bytes):
        assert line, AlignmentError("Empty line")
        assert len(line := decode_line(line).split('\t')) >= 12, AlignmentError(f"Line has < 12 columns: {line}")
        self = cls(
            query_name=line[0], query_length=int(line[1]), query_start=int(line[2]), query_end=int(line[3]),
            strand=line[4], target_name=line[5], target_length=int(line[6]), target_start=int(line[7]),
            target_end=int(line[8]), matching_bases=int(line[9]), num_bases=int(line[10])
        )
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases
        self.percent_query_coverage = 100.0 * (self.query_end - self.query_start) / self.query_length
        self.percent_target_coverage = 100.0 * (self.target_end - self.target_start) / self.target_length
        return self

    def __repr__(self):
        return (f'{self.query_name}:{self.query_start}-{self.query_end} '
                f'{self.target_name}:{self.target_start}-{self.target_end}')

    def __len__(self):
        return self.num_bases

    def target_overlap(self, other: 'Alignment', allowed_overlap: int = 0) -> int:
        """
        Tests whether this alignment overlaps with the other alignment in the query sequence. A bit
        of overlap can be allowed using the allowed_overlap parameter.
        """
        if self is other:
            return 0
        if self.target_name != other.target_name:
            return 0
        # With Minimap2, query_start is always < query_end, so we can skip the sort
        return range_overlap((self.target_start - allowed_overlap, self.target_end + allowed_overlap),
                             (other.target_start, other.target_end), skip_sort=True)

    def query_overlap(self, other: 'Alignment', allowed_overlap: int = 0) -> int:
        """
        Tests whether this alignment overlaps with the other alignment in the target sequence. A bit
        of overlap can be allowed using the allowed_overlap parameter.
        """
        if self is other:
            return 0
        if self.query_name != other.query_name:
            return 0
        # With Minimap2, target_start is always < target_end, so we can skip the sort
        return range_overlap((self.query_start - allowed_overlap, self.query_end + allowed_overlap),
                             (other.query_start, other.query_end), skip_sort=True)

    def conflicts(self, other: 'Alignment', overlap_fraction: float = 0.5):
        """
        Returns whether this hit conflicts with the other hit on the target sequence.
        A conflict is defined as the hits overlapping by 50% or more of the shortest hit's length.
        A hit is not considered to conflict with itself.
        """
        if self is other:
            return False
        if self.target_name != other.target_name:
            return False
        return self.target_overlap(other) / min(self.num_bases, other.num_bases) > overlap_fraction


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=__description__, formatter_class=argparse.RawTextHelpFormatter,
        add_help=False, usage='%(prog)s <reference> <assembly> [options]',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')

    positionals = parser.add_argument_group(bold('Input'), "Breaks are relative to the reference")
    positionals.add_argument('reference', metavar="<gbk>", default=sys.stdin, nargs="?",
                             help='Path to genbank file (default: stdin)', type=argparse.FileType('rt'))
    positionals.add_argument('assembly', help='Path to the (broken) assembly', metavar="<fasta>",
                             type=lambda x: check_file(x))

    options = parser.add_argument_group(bold('Options'))
    options.add_argument('--preset', help='Minimap2 preset (default: %(default)s)', metavar="",
                         default='asm20', choices=['asm5', 'asm10', 'asm20', 'map-ont', 'ava-ont', 'splice'])
    options.add_argument('-m', '--merge', metavar="", default=100, type=int,
                         help='Range merge tolerance (default: %(default)s)')
    options.add_argument('-f', '--feature', metavar="", default="CDS", type=str,
                         help='Feature type to report (default: %(default)s)')
    options.add_argument('-t', '--threads', help='Number of threads to use (default: %(default)s)',
                         metavar="", default=(p := os.cpu_count()), type=lambda x: min(x, p))
    options.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    options.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    return parser.parse_args(a)


def get_feature_name(feature: 'SeqFeature') -> str:
    if "gene" in feature.qualifiers:
        return feature.qualifiers["gene"][0]
    elif "locus_tag" in feature.qualifiers:
        return feature.qualifiers["locus_tag"][0]
    elif feature.id is not None:
        return feature.id
    else:
        return "Unknown feature"


def main():
    args = parse_args(sys.argv[1:])
    check_programs(['minimap2'])

    for r in SeqIO.parse(args.reference, 'genbank'):
        # Align the bad assembly to the reference sequence
        alignments = cull_all_conflicting_alignments(
            minimap2(args.assembly, f'>{r.id}\n{r.seq}\n', preset=args.preset, threads=args.threads))
        if not alignments:
            warning(f"No alignments found between {args.assembly} and {r.id}")
            continue

        for a in alignments:
            for f in r.features:
                if f.type == args.feature:
                    # Check if feature is within the alignment, if so remove it from features
                    if a.target_start <= f.location.start <= f.location.end <= a.target_end:
                        r.features.remove(f)
                    # Check if feature overlaps with the alignment, if so print it
                    elif range_overlap((a.target_start, a.target_end), (f.location.start, f.location.end)) > 0:
                        a_string = f"{a.target_name}\t{a.target_start}\t{a.target_end}\t{a.query_name}\t{a.query_start}\t{a.query_end}"
                        f_name = get_feature_name(f)
                        f_description = f.qualifiers['product'][0] if 'product' in f.qualifiers else ''
                        f_string = f"{f.type}\t{f.location.start}\t{f.location.end}\t{f_name}\t{f_description}"
                        print(f"{a_string}\t{f_string}")
                        r.features.remove(f)

        # remaining features are those that are not in the alignment
        for f in r.features:
            if f.type == args.feature:
                a_string = f"{r.id}\t\t\t\t\t"
                f_name = get_feature_name(f)
                f_description = f.qualifiers['product'][0] if 'product' in f.qualifiers else ''
                f_string = f"{f.type}\t{f.location.start}\t{f.location.end}\t{f_name}\t{f_description}"
                print(f"{a_string}\t{f_string}")


if __name__ == '__main__':
    main()
