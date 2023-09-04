#!/usr/bin/env python
# -*- coding: utf-8 -*

__author__ = 'Tom Stanton'
__title__ = 'Find TFBS'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Find enriched motifs in upstream regions of specified genes'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

import sys
import os
import datetime
from pathlib import Path
from subprocess import Popen, PIPE, call
import argparse
import re
from typing import Generator
from random import randint

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

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


def check_file(file: str | Path) -> Path:
    file = Path(file) if isinstance(file, str) else file
    if not file.is_file():
        quit_with_error(f'{file.name} is not a file')
    elif file.stat().st_size == 0:
        quit_with_error(f'{file.name} is empty')
    else:
        return file.absolute()


def check_dir(path: str | Path) -> Path:
    path = Path(path) if isinstance(path, str) else path
    try:
        path.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        quit_with_error(f"Error creating directory {path}: {e}")
    if not os.access(path, os.R_OK):
        quit_with_error(f"Directory {path} is not readable")
    elif not os.access(path, os.W_OK):
        quit_with_error(f"Directory {path} is not writable")
    return path


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=bold(__description__), formatter_class=argparse.RawTextHelpFormatter,
        add_help=False, usage='%(prog)s <genbank1> <genbank2> ... <genbankN> <gene> [options]',
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')

    positionals = parser.add_argument_group(bold('Input'))
    positionals.add_argument('genbank', nargs="+", help='Genbank files', type=Genome.from_path)
    positionals.add_argument('gene', type=re.compile,
                             help="Regex pattern to match gene feature for TFBS enrichment")
    enrichment_options = parser.add_argument_group(bold('Enrichment options'))
    enrichment_options.add_argument('-u', '--upstream', metavar="", default=300, type=int,
                                    help='Length of upstream region to enrich TFBS (default: %(default)s)')
    enrichment_options.add_argument('-f', '--feature_type', metavar="", default="CDS", type=str,
                                    help='Feature type to search (default: %(default)s)')
    enrichment_options.add_argument('--allowed_genes', metavar="", default=1, type=int,
                                    help='Number of genes per genome for TFBS enrichment (default: %(default)s)')
    enrichment_options.add_argument('--allow_collision', action='store_true',
                                    help='Allow enrichment sequences to overlap with other genes (default: False)')
    enrichment_options.add_argument('--allow_clip', action='store_true',
                                    help='Allow enrichment sequences to be clipped by contig edges (default: False)')

    streme_options = parser.add_argument_group(bold('Streme options'))
    streme_options.add_argument('-p', '--pvalue', metavar="", default=0.05, type=float,
                                help='Streme p-value threshold (default: %(default)s)')
    streme_options.add_argument('-o', '--outdir', metavar="", default="motifs", type=check_dir,
                                help='Output directory for motifs (default: %(default)s)')

    options = parser.add_argument_group(bold('Other options'))
    options.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    options.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                         help='Show program version and exit')

    if len(a) < 2:
        parser.print_help(file=sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


class Genome:
    def __init__(self, path: Path | None = None, enrichment_sequences: list[SeqRecord] | None = None,
                 name: str | None = None, genes_for_enrichment: list[SeqFeature] | None = None):
        self.path = path or Path()
        self.name = name or self.path.stem
        self.enrichment_sequences = enrichment_sequences or []
        self.genes_for_enrichment = genes_for_enrichment or []

    @classmethod
    def from_path(cls, path: str):
        path = check_file(path)
        return cls(path)

    def __iter__(self):
        return self.parse()

    def __repr__(self):
        return self.name

    def parse(self) -> Generator[SeqRecord, None, None]:
        try:
            for record in SeqIO.parse(self.path, 'genbank'):
                yield record
        except Exception as e:
            quit_with_error(f"Error parsing {self.path}: {e}")


def extract_enrichment_sequences(genome: Genome, gene: re.Pattern, feature_type: str, length: int, allowed_genes: int,
                                 allow_collision: bool, allow_clip: bool) -> str:
    log(f"Extracting enrichment sequences from {genome.name}")
    seqs = []
    for contig in genome:
        for f in contig.features:
            if f.type == feature_type:
                if gene.search(str(f.qualifiers)):
                    if f.location.strand == 1:
                        start, end = max(0, f.location.start - length), f.location.start
                    else:
                        start, end = f.location.end, min(len(contig), f.location.end + length)
                    r = contig[start:end] if f.location.strand == 1 else contig[start:end].reverse_complement()
                    r.id = f"{genome.name}_{contig.id}:{start}-{end}_{f.location.strand}"
                    if not allow_clip and len(r) < length:
                        warning(f'Clipped sequence excluded:\n{r}')
                        continue
                    if not allow_collision and len(r.features) > 0:
                        warning(f'Overlapping sequence excluded:\n{r}')
                        continue
                    seqs.append(r)
    if not seqs:
        warning(f"No enrichment sequences found in {genome.name}")
        return ''
    if len(seqs) > allowed_genes:
        warning(f"More than {allowed_genes} genes found in {genome.name}, selecting first {allowed_genes}")
        seqs = seqs[:allowed_genes]
    return ''.join(f">{r.id}\n{r.seq}\n" for r in seqs)


if __name__ == '__main__':
    check_programs(['streme'])
    args = parse_args(sys.argv[1:])

    enrichment_sequences = ''.join(
        extract_enrichment_sequences(g, args.gene, args.feature_type, args.upstream, args.allowed_genes,
                                     args.allow_collision, args.allow_clip) for g in args.genbank)
    if not enrichment_sequences:
        quit_with_error("No enrichment sequences found")

    enrichment_fasta = args.outdir / "enrichment_sequences.fasta"
    enrichment_fasta.write_text(enrichment_sequences)
    seed = randint(0, 100000)
    log(f"Running streme with seed {seed}")
    out, err = Popen(f"streme --seed {seed} --thresh {args.pvalue} --p {enrichment_fasta} --oc {args.outdir}".split(), stdout=PIPE, stderr=PIPE).communicate()
    url = check_file(args.outdir / "streme.html")
    try:  # should work on Windows
        os.startfile(url.as_uri())
    except AttributeError:
        try:  # should work on MacOS and most linux versions
            call(['open', url.as_uri()])
        except Exception as e:
            warning('Could not open html, try opening manually')
    log("Done!")
