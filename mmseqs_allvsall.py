"""
mmseqs_allvsall
-----------------------------------------------------------------------------------------
Function to run an all-vs-all alignment using mmseqs2 (Steinegger and Soding).
-----------------------------------------------------------------------------------------
Requires:
- Python >=3.9
- mmseqs2
- biopython
-----------------------------------------------------------------------------------------
Tom Stanton, 2023
"""

from __future__ import annotations

import os
import sys
from time import time
from subprocess import Popen, PIPE
from pathlib import Path
from argparse import ArgumentParser, RawTextHelpFormatter
from tempfile import TemporaryDirectory

from Bio import SeqIO, SeqRecord


def main():
    check_python_version()
    args = parse_args(sys.argv[1:])
    check_programs(['mmseqs'], verbose=args.verbose)

    start = time()

    # Create a temporary directory to store the concatenated fasta file and the mmseqs2 database
    with TemporaryDirectory() as tmpdir:
        fasta_file = Path(tmpdir) / "concatenated_proteins.faa"
        with open(fasta_file, "wt") as f:
            for file in args.input:
                if file.suffix in [".gbk", ".gbff", ".gb"]:
                    f.write(
                        parse_genbank(file, translate=args.translate, table=args.table, to_stop=args.to_stop,
                                      stop_symbol=args.stop_symbol)
                    )
                elif file.suffix in [".faa", ".fasta", ".fa"]:
                    f.write(parse_fasta(file, table=args.table, to_stop=args.to_stop, stop_symbol=args.stop_symbol))
                else:
                    warning(f"Unknown file type for {file}")

        # Check the concatenated fasta file is not empty
        check_file(fasta_file)

        # Create the database, files are checked inside function
        db_prefix, result_prefix = Path(tmpdir) / "allvsall_db", Path(tmpdir) / "allvsall_result"
        createdb(fasta_file, db_prefix, verbose=args.verbose)

        # Run the alignment
        m8_file = align(db_prefix, db_prefix, result_prefix, threads=args.threads, verbose=args.verbose,
                        extra_args=args.extra_args)

        # Move the result file to the desired location
        m8_file.rename(args.outfile)

        # Optionally keep the temporary directory
        if args.keep:
            Path(tmpdir).rename(Path.cwd() / f"mmseqs_allvsall_{time()::0.0f}")

    # Print the time taken
    log(f"Result file written to {args.outfile} ({time() - start:.2f} seconds)", end="\n\n")


DESCRIPTION = """
Run mmseqs2 all-vs-all alignment
-----------------------------------------------------------------------------------------
    Takes fasta or genbank files containing CDS/protein sequences and 
    runs an all-vs-all alignment using mmseqs2.
-----------------------------------------------------------------------------------------
"""


def parse_args(args):
    parser = ArgumentParser(description=DESCRIPTION, formatter_class=RawTextHelpFormatter,
                            epilog="Note: mmseqs2 must be installed and in your PATH",
                            usage="python3 mmseqs_allvsall.py file.{fasta,gbk} [file.{fasta,gbk}... ] out.m8")
    parser.add_argument("input", nargs="+", type=lambda x: check_file(Path(x)),
                        help="Fasta or Genbank file(s) containing the CDS/Protein sequences to align.\n"
                             "If a DNA fasta file is provided, it will be translated on the fly.\n"
                             "If a Genbank file is provided, the translation in the Genbank features will be used\n"
                             "Unless --translate is used, in which case the DNA sequence will be translated.\n")
    parser.add_argument("outfile",
                        help="Path for the M8-formatted tabular result file", type=lambda x: Path(x).absolute())
    parser.add_argument("--extra_args",
                        help="Extra arguments to pass to mmseqs2", default='', metavar='[None]')
    parser.add_argument("--translate", action="store_true", default=False,
                        help="DNA features will be extracted and translated with BioPython.\n"
                             "If False, the translation in the Genbank features will be used.")
    parser.add_argument("--table", default=11, help="Translation table for BioPython translate", metavar='[11]')
    parser.add_argument("--to_stop", action="store_true", default=False,
                        help="Stop translation at first stop codon for BioPython translate")
    parser.add_argument("--stop_symbol", type=str, default="*", metavar="[*]",
                        help="Stop symbol for translation")
    parser.add_argument("--keep", action="store_true", default=False,
                        help="Move temporary directory to CWD after run")
    parser.add_argument("--feature_type", type=str, default="CDS", metavar="[CDS]",
                        help="Genbank feature type to extract")
    parser.add_argument("-t", "--threads", default=(cpus := os.cpu_count()), metavar=f"[{cpus}]",
                        type=lambda x: min(int(x), cpus), help="Number of threads to use")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    return parser.parse_args(args)


def check_file(path: Path):
    if not path.is_file():
        quit_with_error(f"File {path} does not exist")
    if not path.stat().st_size:
        quit_with_error(f"File {path} is empty")
    return path.absolute()


def log(message: str, flush: bool = True, out=sys.stderr, end='\n'):
    print(message, flush=flush, file=out, end=end)


def warning(message: str, **kwargs):
    message = f"\033[1;33m{message}\033[0m"  # Wrap message in bold yellow text
    log(f"WARNING: {message}", **kwargs)


def error(message: str, **kwargs):
    message = f"\033[1;31m{message}\033[0m"  # Wrap message in bold red text
    log(f"ERROR: {message}", **kwargs)


def quit_with_error(message: str, **kwargs):
    error(message, **kwargs)
    sys.exit(1)


def check_python_version(major=3, minor=9):
    if sys.version_info.major < major or sys.version_info.minor < minor:
        quit_with_error("Python 3.9 or greater is required")


def align(query_prefix: Path, target_prefix: Path, result_prefix: Path, extra_args: str = '', threads: int = 1,
          verbose: bool = False):
    if query_prefix == target_prefix:
        if verbose:
            log(f"Creating fake prefix for all-vs-all alignment")
        create_fake_prefix(target_prefix, result_prefix)

    tsv_path = result_prefix.with_suffix(".m8")
    alignment_prefix = Path(str(result_prefix) + '_aln')

    run_command(
        f'mmseqs align {query_prefix} {target_prefix} {result_prefix} {alignment_prefix} --threads {threads} {extra_args}'.strip(),
        verbose=verbose)

    run_command(f'mmseqs convertalis {query_prefix} {target_prefix} {alignment_prefix} {tsv_path}', verbose=verbose)

    return check_file(tsv_path)


def create_fake_prefix(target_prefix: Path, allvsall_prefix: Path):
    """create a fake prefiltering for all-vs-all alignments"""

    # create link to data file which contains a list of all targets that should be aligned
    target_db_index = check_file(target_prefix.with_suffix(".index"))
    allvsall_prefix.symlink_to(target_db_index)
    allvsall_index = allvsall_prefix.with_suffix(".index")
    index_size = target_db_index.stat().st_size

    # create new index repeatedly pointing to same entry
    with open(target_db_index, "rt") as f, open(allvsall_index, "wt") as g:
        g.write("\n".join(line.split("\t")[0] + "\t0\t" + str(index_size) for line in f.readlines()) + "\n")

    # create dbtype (7)
    allvsall_dbtype = allvsall_prefix.with_suffix(".dbtype")
    allvsall_dbtype.write_text(f"{chr(7)}{chr(0)}{chr(0)}{chr(0)}")
    return check_file(allvsall_dbtype)


def createdb(fasta: Path | list[Path], prefix: Path, extra_args: str = "", verbose: bool = False):
    """
    Create a MMseqs2 db
    """
    command = f"mmseqs createdb {extra_args}".strip()
    command += f" {' '.join([str(i) for i in fasta])} {prefix}" if isinstance(fasta, list) else f" {fasta} {prefix}"
    run_command(command, verbose=verbose)
    return find_files_with_suffixes(prefix, ['.dbtype', '.index', '.lookup', '.source'])


def find_files_with_suffixes(prefix: Path, suffixes: list[str], min_size: int = 1) -> list[Path]:
    """
    Find files with given suffixes for the given prefix, useful for finding blast database or bwa index files
    :param suffixes: List of existing Paths with suffixes
    :param prefix: Prefix as Path
    :param min_size: Minimum size of file in bytes
    :return: List[Path]
    """
    return [p for i in suffixes if
            (p := prefix.with_suffix(prefix.suffix + i)).exists() and p.stat().st_size >= min_size]


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
            quit_with_error(f'{program} not found in PATH')


def run_command(command: str, pipe: bytes = b'', shell: bool = False, cmd_split: str = ' ', verbose: bool = False,
                string_out: bool = False):
    """Run a command and return the output"""
    start = time()
    if verbose:
        log(f'Running Command\n\tPipe: {"True" if pipe else "False"}\n'
            f'\tShell: {"True" if shell else "False"}\n\tCommand: {command}')
        log(f'\r\tRunning...', flush=False, end=' ')
    if pipe:
        with Popen(command.split(cmd_split), shell=shell, stdout=PIPE, stderr=PIPE, stdin=PIPE) as child:
            out, err = child.communicate(input=pipe)
    else:
        with Popen(command.split(cmd_split), shell=shell, stdout=PIPE, stderr=PIPE) as child:
            out, err = child.communicate()

    if verbose:
        log(f'Done! ({time() - start:.2f}s)')
        if err and not out:
            warning(err.decode())

    return out.decode() if string_out else out


def parse_fasta(fasta: Path, to_stop: bool = True, table: int = 11, stop_symbol: str = '*'):
    """Parse a fasta file and write to a file handle"""
    protein_records = ''
    try:
        for record in SeqIO.parse(fasta, 'fasta'):
            if set(record.seq).issubset(set('ACGTN')):  # Check if seq is DNA
                protein_records += f'>{record.id}\n' \
                                   f'{record.seq.translate(to_stop=to_stop, table=table, stop_symbol=stop_symbol)}\n'
            else:
                protein_records += f'>{record.id}\n{record.seq}\n'
    except ValueError:
        warning(f'Failed to parse {fasta}')
    return protein_records


def parse_genbank(genbank: Path, translate: bool = True, to_stop: bool = True, table: int = 11, stop_symbol: str = '*'):
    """Parse a genbank file and write to a file handle"""
    protein_records = ''
    try:
        protein_records += ''.join(parse_gb_record(r, "CDS", translate, to_stop, table, stop_symbol) for r in SeqIO.parse(genbank, 'genbank'))
    except ValueError:
        warning(f'Failed to parse {genbank}')
    return protein_records


def parse_gb_record(record: SeqRecord, feature_type: str = "CDS", translate: bool = True, to_stop: bool = True,
                    table: int = 11, stop_symbol: str = '*'):
    protein_records = ''
    for feature in record.features:
        if feature.type == feature_type:
            identifier = feature.qualifiers['locus_tag'][0]
            if translate:
                protein_records += f">{identifier}\n" \
                                   f"{feature.extract(record.seq).translate(to_stop=to_stop, table=table, stop_symbol=stop_symbol)}\n"
            else:
                if 'translation' in feature.qualifiers:
                    protein_records += f">{identifier}\n" \
                                       f"{feature.qualifiers['translation'][0]}\n"
                else:
                    warning(f'No translation for {identifier}')
                    continue
    if not protein_records:
        warning(f'No {feature_type} found in {record.id}')
    return protein_records


if __name__ == '__main__':
    main()
