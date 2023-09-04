#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Tom Stanton'
__title__ = 'Gluggins'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Parse Gubbins output as recombinant blocks of features'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import sys
import datetime
from pathlib import Path
from typing import Generator
import argparse

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

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


def check_file(file: str | Path) -> Path:
    file = Path(file) if isinstance(file, str) else file
    if not file.is_file():
        quit_with_error(f'{file.name} is not a file')
    elif file.stat().st_size == 0:
        quit_with_error(f'{file.name} is empty')
    else:
        return file.absolute()


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
    return max(0, min(end1, end2) - max(start1, start2))


def blocks_from_events(events: 'GubbinsGff', distance: int = 0, min_nll: float = 0) -> Generator['Block', None, None]:
    """
    Piles up regions of overlapping events into blocks
    """
    log(f"Generating blocks from overlapping events with a distance of {distance} bp and a minimum NLL of {min_nll}...")
    events_in_block = []
    for event in sorted(events, key=lambda i: i.start):
        if event.neg_log_likelihood >= min_nll:
            if len(events_in_block) == 0 or event.start <= events_in_block[-1].end + distance:
                events_in_block.append(event)
            else:
                yield Block(events=events_in_block, start=events_in_block[0].start, end=events_in_block[-1].end)
                events_in_block = []


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=bold(__description__), formatter_class=argparse.RawTextHelpFormatter,
        add_help=False, usage='%(prog)s <reference> <gff> [options] > out.{tab,gbk,fa}',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')

    opts = parser.add_argument_group(bold('Input'),
                                     """
Gluggins assumes first record in reference corresponds to the sequence used for Gubbins.
An error will be raised if the length of the reference does not match the sequence length in the Gubbins gff.
""")
    opts.add_argument('reference', help='Genbank file or - for stdin', type=argparse.FileType('rt'))
    opts.add_argument('gff', help='Path to Gubbins gff', type=GubbinsGff.from_path)

    opts = parser.add_argument_group(bold('Output options'),
                                     """
Gluggins can output blocks of features in tab-delimited, genbank or fasta format to stdout.
The fasta format is a multi-fasta with each block as a separate record.
The genbank format is a single genbank with each block as a separate record adding events as features.
The default output is a tab-delimited file with the following columns:

    1. Reference sequence ID
    2. Block ID
    3. Block start position
    4. Block end position
    5. Feature type
    6. Feature start position
    7. Feature end position
    8. Feature strand
    9. Feature name
    10. Feature locus tag
    11. Feature description
    12. Events
    
    The events column is a semi-colon-delimited list of all events overlapping each feature in the block.
    The events are space-delimited with the following columns:
    
        1. Start position
        2. End position
        3. Node ID
        4. SNP count
        5. A comma-delimited list of taxa
""")

    opts.add_argument('-f', '--format', choices=['tab', 'genbank', 'fasta'], default='tab', metavar="",
                      help='Output format [%(choices)s] (default: %(default)s)')
    opts.add_argument('-t', '--feature_type', metavar="", default='CDS',
                      help='Feature type to extract from reference (default: %(default)s)')
    opts.add_argument('-s', '--suppress_header', action='store_true',
                      help='Suppress header in tab-delimited output')

    opts = parser.add_argument_group(bold('Block options'), '')
    opts.add_argument('-d', '--distance', metavar="", default=1000, type=int,
                      help='Maximum distance between events in the same block (default: %(default)s)')
    opts.add_argument('-n', '--min_nll', metavar="", default=0, type=float,
                      help='Minimum block negative log likelihood (default: %(default)s)')
    opts.add_argument('-e', '--min_events', metavar="", default=1, type=int,
                      help='Minimum number of events per block (default: %(default)s)')
    opts.add_argument('-l', '--min_len', metavar="", default=100, type=int,
                      help='Minimum length per block (default: %(default)s)')
    opts.add_argument('-m', '--min_features', metavar="", default=1, type=int,
                      help='Minimum number of features per block (default: %(default)s)')

    other_options = parser.add_argument_group(bold('Other options'), '')
    other_options.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    other_options.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                               help='Show program version and exit')

    if len(a) < 2:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


# Classes ----------------------------------------------------------------------
class GubbinsGff:
    def __init__(self, path: Path | None = None, name: str | None = None, sequence_length: int | None = None,
                 events: list['Event'] | None = None):
        self.path = path or Path()
        self.name = name or self.path.stem
        self.sequence_length = sequence_length or 0
        self.events = events or []

    @classmethod
    def from_path(cls, path: str):
        path = check_file(path)
        self = cls(path=path, name=path.stem)
        with open(path, 'rt') as fh:
            line_1 = fh.readline()
            if not line_1.startswith('##gff-version 3'):
                quit_with_error(f'{self} does not start with ##gff-version 3')
            line_2 = fh.readline()
            if not line_2.startswith('##sequence-region'):
                quit_with_error(f'{self} does not have a sequence region')
            self.sequence_length = int(line_2.split()[3])
            for line in fh.readlines():
                if not line.startswith('#'):
                    try:
                        self.events.append(Event.from_gff(line))
                    except EventError as e:
                        warning(f"Skipping line: {e}")
        if not self.events:
            quit_with_error(f"No events found in {self}")
        return self

    def __iter__(self):
        return iter(self.events)

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.events)


class EventError(Exception):
    pass


class Event:
    def __init__(
            self, seqid: str | None = None, source: str | None = None, type_: str | None = None,
            start: int | None = None, end: int | None = None, score: float | None = 0, strand: str | None = None,
            phase: float | None = 0, node: str | None = None, net_log_likelihood: float | None = None,
            taxa: list[str] | None = None, snp_count: int | None = None, features: list[SeqFeature] | None = None,
            record: SeqRecord | None = None
    ):
        self.seqid = seqid or ''
        self.source = source or ''
        self.type = type_ or ''
        self.start = start or 0
        self.end = end or 0
        self.score = score or 0
        self.strand = strand or ''
        self.phase = phase or 0
        self.node = node or ''
        self.neg_log_likelihood = net_log_likelihood or 0
        self.taxa = taxa or []
        self.snp_count = snp_count or 0
        self.features = features or []
        self.record = record or None

    @classmethod
    def from_gff(cls, line: str | bytes):
        if not line:
            raise EventError("Empty line passed to Event")
        line = line.decode().strip() if isinstance(line, bytes) else line.strip()
        if line.startswith('#'):
            raise EventError("Comment line passed to Event")
        if len(line := line.split('\t')) < 9:
            raise EventError(f"Gff line has less than 9 columns: {line}")
        if line[1] != 'GUBBINS':
            raise EventError(f"Unexpected source in GFF: {line[1]}")

        self = cls(
            seqid=line[0], source=line[1], type_=line[2], end=int(line[4]), score=float(line[5]), strand=line[6],
            phase=float(line[7]), start=int(line[3]) - 1)  # Convert to 0-based

        for attr in line[8].split(';'):
            if attr:
                if len(i := attr.split('=', 1)) == 2:
                    if i[0] == 'node':
                        self.node = i[1].replace('"', '')
                    elif i[0] == 'neg_log_likelihood':
                        self.neg_log_likelihood = float(i[1].replace('"', ''))
                    elif i[0] == 'taxa':
                        self.taxa = [i for i in i[1].replace('"', '').split(' ') if i]
                    elif i[0] == 'snp_count':
                        self.snp_count = int(i[1].replace('"', ''))
                    else:
                        raise EventError(f"Unexpected attribute in GFF: {i}")
                else:
                    raise EventError(f"Could not parse attribute in GFF: {attr}")
        return self

    def __len__(self):
        return self.end - self.start

    def __repr__(self):
        return f"{self.start}-{self.end}"

    def as_line(self):
        return f"{self.start} {self.end} {self.node} {self.snp_count} {','.join(self.taxa)}"

    def as_seqfeature(self) -> SeqFeature:
        return SeqFeature(
            type="misc_recomb", location=SimpleLocation(self.start, self.end, 1 if self.strand == "+" else -1),
            id=f"gubbins_event_{self}",
            qualifiers={
                "inference": "GUBBINS",
                "note": f"node={self.node};neg_log_likelihood={self.neg_log_likelihood};taxa={','.join(self.taxa)};"
                        f"snp_count={self.snp_count}"
            }
        )


class Block:
    """
    Class to represent a block of overlapping events and features
    """
    def __init__(self, start: int | None = None, end: int | None = None,
                 events: list[Event] | None = None, features: list['Feature'] | None = None):
        self.start = start or 0
        self.end = end or 0
        self.events = events or []
        self.features = features or []

    def __repr__(self):
        return ','.join(str(i) for i in self.features) if self.features else ""

    def __len__(self):
        return self.end - self.start

    def add_features(self, record: SeqRecord, feature_type: str):
        """
        Add a SeqRecord to the block and extract features that overlap the block.
        Note, simply slicing the record changes the relative location of the features, so we need to iterate over them
        """
        self.features += [Feature.from_seqfeature(f, self) for f in record.features if f.type == feature_type and
                          range_overlap((self.start, self.end), (f.location.start, f.location.end)) > 0]
        # Extend the block to include the features
        self.start = min(self.start, min(f.start for f in self.features))
        self.end = max(self.end, max(f.end for f in self.features))

    def as_seqrecord(self, record: SeqRecord) -> SeqRecord:
        """
        Return a SeqRecord containing the block features and events
        """
        record.features += [e.as_seqfeature() for e in self.events]  # Add the events as features
        return record[self.start:self.end]  # Slice the record to the block


class Feature:
    """
    Class to represent a feature (default: CDS) that overlaps at least one event
    """

    def __init__(self, feature: SeqFeature | None = None, events: list[Event] | None = None, block: Block | None = None,
                 name: str | None = None, locus_tag: str | None = None, description: str | None = None,
                 start: int | None = None, end: int | None = None, strand: str | None = None):
        self.feature = feature or None
        self.block = block or None
        self.name = name or ''
        self.start = start or 0
        self.end = end or 0
        self.locus_tag = locus_tag or ''
        self.description = description or ''
        self.strand = strand or ''
        self.events = events or []

    @classmethod
    def from_seqfeature(cls, feature: SeqFeature, block: Block):
        return cls(
            feature=feature, block=block, name=feature.qualifiers.get('gene', [''])[0],
            locus_tag=feature.qualifiers.get('locus_tag', [''])[0],
            description=feature.qualifiers.get('product', [''])[0],
            strand="+" if feature.strand == 1 else "-",
            start=(start := int(feature.location.start)), end=(end := int(feature.location.end)),
            events=[e for e in block.events if range_overlap((start, end), (e.start, e.end)) > 0]
        )

    def __len__(self):
        return self.end - self.start

    def __repr__(self):
        return self.name if self.name else self.locus_tag if self.locus_tag else self.description

    def as_line(self):
        """
        Note, converting to 1-based
        """
        return (f"{self.start + 1}\t{self.end}\t{self.strand}\t{self.name}\t{self.locus_tag}\t{self.description}\t"
                f"{';'.join(i.as_line() for i in self.events)}")


if __name__ == '__main__':
    # def main():
    args = parse_args(sys.argv[1:])
    # args = parse_args(["SL_10_KP_NORM_BLD_13379_reference.gbff", 'SL_10_gubbins.recombination_predictions.gff'])
    args.reference = SeqIO.read(args.reference, 'genbank')
    if args.gff.sequence_length != len(args.reference):
        quit_with_error(f"{args.reference.id} is {len(args.reference)} bp, is {args.gff.sequence_length} bp")

    blocks = []
    for b in blocks_from_events(args.gff, args.distance, args.min_nll):
        if len(b.events) >= args.min_events and len(b) >= args.min_len:
            b.add_features(args.reference, args.feature_type)
            if len(b.features) >= args.min_features:
                blocks.append(b)
    if not blocks:
        quit_with_error("No blocks found passing filters")

    log(f"{len(blocks)} blocks passing filters: min_events={args.min_events}, min_len={args.min_len}, "
        f"min_features={args.min_features}")

    if args.format == 'fasta':
        sys.stdout.write(
            ''.join(b.as_seqrecord(args.reference).format('fasta') for b in blocks)
        )

    elif args.format == 'genbank':
        sys.stdout.write(
            ''.join(b.as_seqrecord(args.reference).format('genbank') for b in blocks)
        )

    else:
        if not args.suppress_header:
            sys.stdout.write(
                f"reference\tblock_id\tblock_start\tblock_end\tfeature_type\tfeature_start\tfeature_end\t"
                f"feature_strand\tfeature_name\tfeature_locus_tag\tfeature_description\tevents\n")

            sys.stdout.write(
                ''.join(  # Note, converting to 1-based
                    ''.join(f"{args.reference.id}\tBlock_{n}\t{b.start + 1}\t{b.end}\t{args.feature_type}\t"
                            f"{f.as_line()}\n" for f in b.features) for n, b in enumerate(blocks, 1)
                )
            )

    log(f"Done!")
