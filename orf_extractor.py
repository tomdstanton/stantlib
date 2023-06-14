#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Tom Stanton'
__title__ = 'ORF Extractor'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Quickly extract ORFs in bacterial genomes'
__long_description__ = 'Annotate IS elements and interrupted genes in bacterial genomes'
__license__ = 'gpl-3.0'
__version__ = '0.0.1'

# ................ Python Imports ................ #
import logging
import os
from pathlib import Path
import sys
import gzip
import hashlib
from typing import Union, List
from textwrap import wrap
import multiprocessing.pool
from operator import attrgetter
from argparse import ArgumentParser, RawTextHelpFormatter, SUPPRESS
# ................ Dependency Imports ................ #
import pyrodigal

# ................ Set global variables ................ #
LINE_WIDTH = 70  # Nice general terminal width
LOGGER = logging.getLogger(__name__)  # https://stackoverflow.com/a/45065655/10771025
ASCII = fr'''
{__title__}
{__description__}
'''


# ................ Classes ................ #
class Assembly:
    def __init__(self, path: Path, arguments: 'argparse.Namespace'):
        self.protein_hits = None
        self.name = path.stem
        self.path = path
        self.extension = path.suffix
        self.genes = {}
        self.args = arguments
        self.orf_finder = pyrodigal.OrfFinder(
            meta=False, closed=False, mask=False, min_gene=self.args.min_gene,
            min_edge_gene=self.args.min_edge_gene, max_overlap=self.args.max_overlap
        )
        self.contigs = {}
        if arguments.locus_tag:
            self.locus_tag = arguments.locus_tag
        else:
            self.locus_tag = generate_locus_tag(path) if arguments.auto_locus_tag == 'prokka' else self.name

    def __repr__(self):
        return self.name

    def __len__(self):
        return self.path.stat().st_size

    def largest_contig(self):
        return max(self.contigs.values(), key=attrgetter('length'))

    def get_genes(self):
        LOGGER.info(f"Finding genes in {self.name}")
        # https://pyrodigal.readthedocs.io/en/stable/api/orf_finder.html#pyrodigal.OrfFinder.train
        if self.args.train == 'concatenated':  # Train on concatenated contigs
            self.orf_finder.train(
                'TTAATTAATTAA'.join([i.seq for i in self.contigs.values()]), translation_table=self.args.translation_table)
        else:
            self.orf_finder.train(
                self.largest_contig().seq, translation_table=self.args.translation_table)  # Train on largest contig
        with multiprocessing.pool.ThreadPool(self.args.threads) as pool:
            predicted_genes = pool.map(self.orf_finder.find_genes, [contig.seq for contig in self.contigs.values()])
        if predicted_genes:
            gene_counter = 0
            for contig, genes in zip(self.contigs.values(), predicted_genes):
                n_genes_on_contig = len(genes)
                for position, gene in enumerate(genes):
                    contig.genes.append(CDS(gene, contig, position, n_genes_on_contig, gene_counter))
                    gene_counter += 1
            LOGGER.info(f"Found {gene_counter + 1} genes in {self.name}")
        else:
            LOGGER.warning(format_str(f'No genes found in {self.name}', '93'))

    def write_genes(self):
        if self.genes:
            for mol in self.args.molecule:
                if self.args.out:
                    outfile = Path(self.args.out / f'{self.name}_{mol}.fasta')
                    with open(outfile, 'wt') as f:
                        f.write('\n'.join(
                            [i.get_fasta(mol=mol) for i in self.genes.values()]
                        ) + '\n')
                    if outfile.stat().st_size > 0:
                        LOGGER.info(f'Wrote {outfile}')
                    else:
                        LOGGER.warning(format_str(f'{outfile} empty, deleting', '93'))
                        outfile.unlink()
                else:
                    print('\n'.join(
                        [i.get_fasta(mol=mol) for i in self.genes.values()]
                    ))


class AssemblyFasta(Assembly):
    def __init__(self, path, arguments, fasta_records: List['SeqRecord']):
        super().__init__(path, arguments)
        self.contigs = {record.name: Contig(record, record_n + 1, self) for record_n, record in
                        enumerate(fasta_records)}


class SeqRecord:
    def __init__(self, fasta_lines: List[str]):
        self.name, self.description = fasta_lines[0].split(' ', 1) if ' ' in fasta_lines[0] else (fasta_lines[0], '')
        self.sequence = ''.join(fasta_lines[1:])

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.sequence)


class Contig:
    """
    Class for storing CDS predicted by P(y)rodigal in their relative positions on a sequence.
    They are stored in a list and can be accessed by their relative position and the 'CDS.list_lookup' attribute.
    """

    def __init__(self, record, contig_number: int, a: Assembly):
        self.record = record
        self.assembly = a
        self.name = record.name
        self.seq = record.sequence
        self.length = len(self.seq)
        self.genes = []
        self.start = 0  # Zero based
        self.end = self.length  # Zero based
        self.contig_number = contig_number

    def __repr__(self):
        return self.name

    def __len__(self):
        return self.length

    def __eq__(self, other: 'Contig'):
        return self.name == other.name


class CDS:
    """
    Class for CDS predicted by P(y)rodigal for interacting with a single gene within a contig, within an assembly.
    This class is used by many stantlib programs and should be relatively unmodified.
    NOTE: P(y)rodigal uses 1-based coordinates, so the start and end positions are 1-indexed.
    """

    def __init__(self, g: pyrodigal.Gene, c: Contig, position_in_contig, last_position_on_contig, position_in_assembly):
        self.contig = c
        self.assembly = c.assembly
        self.assembly_type = self.assembly.__class__.__name__
        self.list_lookup = position_in_contig  # Zero based
        self.relative_position_on_contig = position_in_contig + 1  # One based
        self.relative_position_on_assembly = position_in_assembly + 1  # One based
        self.gene = g
        self.strand = "+" if g.strand == 1 else "-"
        self.start = g.begin - 1  # Zero based
        self.end = g.end - 1  # Zero based
        self.name = f'{self.assembly.locus_tag}_{str(self.relative_position_on_assembly).zfill(5)}'
        self.assembly.genes[self.name] = self
        self.edge_of_contig = True if self.relative_position_on_contig in [1, last_position_on_contig] else False

    def __repr__(self):
        return f'{self.name} Contig: {self.contig.name} ({len(self.contig.genes)} CDS)'

    def __len__(self):
        return self.gene.end - self.gene.begin

    def __eq__(self, other: 'CDS'):
        return self.name == other.name

    def partial(self):
        return any([self.gene.partial_begin, self.gene.partial_end])

    def only_cds_on_contig(self):
        return len(self.contig.genes) == 1

    def get_fasta(self, mol=Union['protein', 'dna'], seq_wrap=0) -> str:
        """ Returns the fasta-formatted sequence of the gene in the specified molecule"""
        seq = self.gene.translate() if mol == 'protein' else self.gene.sequence()
        return f'>{self.name}\n{fasta_wrap(seq, seq_wrap) if seq_wrap else seq}'


# ............... Functions ............... #
def check_file(parser: 'argparse.ArgumentParser', file: Path) -> Path:
    if not file.is_file():
        parser.error(f'{file} does not exist')
    elif file.stat().st_size == 0:
        parser.error(f'{file} is empty')
    else:
        return file


def check_dir(parser: 'argparse.ArgumentParser', dirpath) -> Path:
    try:
        dirpath.mkdir(parents=True, exist_ok=True)
        return dirpath
    except Exception as e:
        parser.error(str(e))


def determine_file_type(file) -> Union[str, None]:
    try:
        first_byte = file.read(1)
        if first_byte == '>':
            return 'fasta'
        elif first_byte == 'S':
            return 'gfa'
    except Exception as e:
        LOGGER.warning(format_str(str(e), '93'))
        return None


def generate_locus_tag(file: Path, ) -> str:
    # https://stackoverflow.com/a/3431835/10771025
    # https://github.com/tseemann/prokka/blob/c335b7a6539a960ff4d98af6d2c6eb9bc1b7ebde/bin/prokka#L1702
    hash_function = hashlib.md5()
    with open(file, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_function.update(chunk)
    # chr(65) = A
    return ''.join([chr(65 + s) for s in [int(i) for i in hash_function.hexdigest() if i.isdigit()][0:8]])


def check_assembly(path: Path, arguments: 'argparse.Namespace'):
    LOGGER.debug(f'Checking assembly file: {path}')
    try:
        with open(path, 'rt') as f:
            data = f.read()
    except IOError as e:
        LOGGER.warning(format_str(str(e), '93'))
        try:
            with gzip.open(path, 'rt') as f:
                data = f.read()
        except Exception as e:
            LOGGER.warning(format_str(str(e), '93'))
            LOGGER.error(format_str(f'Could not open {path}', '91'))
            return None

    if data[0] == '>':
        try:
            return AssemblyFasta(path, arguments, [SeqRecord(i.splitlines()) for i in data.split('>')[1:] if i])
        except Exception as e:
            LOGGER.warning(format_str(f'Could not create AssemblyFasta object: {e}', '93'))
            return None


def fasta_wrap(fasta_string: str, width: int = 80) -> str:
    return '\n'.join(wrap(fasta_string, width))


def format_str(string: str, format_number: str) -> str:
    """Add colour to string using a format string number"""
    return f"\033[{format_number}m{string}\033[0m"


def quit_with_error(message: str):
    """Displays the given message and ends the program's execution."""
    LOGGER.error(format_str(message, '91'))
    sys.exit(1)


def parse_args(arguments: list) -> 'argparse.Namespace':
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter, description=format_str(ASCII, '36'),
                            usage='python3 interupta.py [options] [assembly.fasta(.gz)]', add_help=False)
    parser.add_argument('input', nargs='+', help=SUPPRESS,
                        type=lambda p: check_file(parser, Path(p).resolve()))

    orf_parser = parser.add_argument_group('ORF options')
    orf_parser.add_argument('--min_gene', type=int, default=90, help="The minimum gene length")
    orf_parser.add_argument('--min_edge_gene', type=int, default=60,
                            help="The minimum edge gene length")
    orf_parser.add_argument('--max_overlap', type=int, default=60,
                            help="The maximum number of nucleotides that can overlap between two genes on the same strand")
    orf_parser.add_argument('--train', choices=['largest', 'concatenated'], default='largest',
                            help="Pyrodigal contig training method")
    orf_parser.add_argument('--translation_table', type=int, default=11, help="Translation table to use")
    orf_parser.add_argument('--auto_locus_tag', choices=['filename', 'prokka'], default='filename',
                            help="Automatic locus tag method")
    orf_parser.add_argument('--locus_tag', type=str, nargs="?", help="Custom locus tag")

    output_parser = parser.add_argument_group('Output options')
    output_parser.add_argument('-o', '--out', type=lambda p: check_dir(parser, Path(p).resolve()),
                               const=Path.cwd(), nargs="?",
                               help="Write hit and interrupted fasta files in specified directory")
    output_parser.add_argument('-m', '--molecule', choices=['protein', 'dna'], default=['protein'], nargs="+",
                               help="Molecule type for hit and interrupted fasta files")

    other_parser = parser.add_argument_group('Other options')
    other_parser.add_argument('-t', '--threads', type=lambda t: min([os.cpu_count(), int(t)]), default=os.cpu_count(),
                              help="Number of threads to use")
    other_parser.add_argument('--version', action='version', version=f'{__name__} v{__version__}', help=SUPPRESS)
    other_parser.add_argument('-v', '--verbosity', action="count", default=0, help="verbosity, increases with -vv")
    other_parser.add_argument('-h', '--help', action='help', help='show this help message and exit')
    return check_args(parser.parse_args(arguments))


def check_args(arguments: 'argparse.Namespace') -> 'argparse.Namespace':
    if arguments.verbosity == 1:
        arguments.verbosity = logging.INFO
    elif arguments.verbosity > 1:
        arguments.verbosity = logging.DEBUG
    else:
        arguments.verbosity = logging.WARNING
    logging.basicConfig(format="[%(asctime)s %(levelname)s %(funcName)s] %(message)s",
                        datefmt='%H:%M:%S', level=arguments.verbosity)
    LOGGER.info(f' {Path(__file__).stem} {__version__} '.center(LINE_WIDTH, '='))
    LOGGER.info(f'Platform: {sys.platform}')
    LOGGER.info(f'Python: {".".join([str(i) for i in sys.version_info[:3]])}')
    LOGGER.info(f'Current working directory: {Path.cwd()}')
    if not arguments.input:
        quit_with_error('No assembly fasta provided')
    return arguments


def progress_bar(current: int, total: int, bar_length: int = 20) -> str:
    """Progress bar for terminal output"""
    percent = float(current) * 100 / total
    arrow = '=' * int(percent / 100 * bar_length - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    col = '96' if percent < 100 else '92'
    return format_str(f'[{arrow}{spaces}] {percent:.1f}%', col)


if __name__ == '__main__':
    # ................ Parse Args ................ #
    args = parse_args(sys.argv[1:])  # Uncomment below for testing
    # args = parse_args(['/Users/tom/Bioinformatics/stantlib/ERR4920398.fasta'])

    # ................ Iterate over assemblies ................ #
    for n, file in enumerate(args.input):
        assembly = check_assembly(file, args)
        if assembly:
            assembly.get_genes()
            assembly.write_genes()
            print(progress_bar(n + 1, len(args.input), LINE_WIDTH), end="\r", flush=True, file=sys.stderr)
        else:
            LOGGER.error(format_str(f"Could not parse input file: {file}", '91'))
    sys.exit(0)
