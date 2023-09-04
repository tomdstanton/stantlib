#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Tom Stanton'
__title__ = 'Get HMMs'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Download Hidden-Markov/Covariance models from public databases with accessions'
__license__ = 'gpl-3.0'
__version__ = '0.0.1b'

from pathlib import Path
import requests
from zlib import decompress, MAX_WBITS
import datetime
import sys
from os import cpu_count
import argparse
from multiprocessing.dummy import Pool

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


def download(url: str) -> bytes:
    hmm_string = b''
    log(f"Downloading {url}")
    with requests.get(url, stream=True) as request:
        if request.status_code == 200:
            for chunk in request.iter_content(chunk_size=2 ** 20):
                hmm_string += chunk
            hmm_string = decompress(hmm_string, wbits=MAX_WBITS | 32) if is_gzipped(hmm_string) else hmm_string
            return hmm_string
        else:
            warning(f"Could not download {url}, status code {request.status_code}")
    return hmm_string


def check_cpus(cpus: int | str) -> int:
    try:
        cpus = int(cpus)
    except ValueError:
        quit_with_error(f"CPUs must be an integer, got {cpus}")
    if cpus < 1:
        quit_with_error(f"CPUs must be > 0, got {cpus}")
    return min(cpus, cpu_count())


def parse_args(a):
    parser = argparse.ArgumentParser(
        description=bold(__description__), add_help=False,
        usage='%(prog)s <accessions> [options]', formatter_class=argparse.RawTextHelpFormatter,
        epilog=f'Author: {__author__}\tEmail: {__author_email__}\tLicense: {__license__}\tVersion: {__version__}')

    opts = parser.add_argument_group(bold('Input'), """
Currently supported databases:
 - PFAM: PF00001
 - RFAM: RF00001
 - CDD: cd00001
 - PANTHER: PTHR00001
 - TIGRFAM/NCBIFam: TIGR00001 or NF00001
""")
    opts.add_argument('accessions', nargs='+', help="HMM/CM accessions to download")

    opts = parser.add_argument_group(bold('Options'))
    opts.add_argument('-o', '--outdir', default=Path.cwd(), type=Path, metavar='',
                      help='Output directory (default: cwd)')
    opts.add_argument('-e', '--extension', default='hmm', type=str, metavar='',
                      help='HMM/CM file extension (default: hmm)')
    opts.add_argument('-t', '--threads', default=cpu_count(), type=check_cpus, metavar='',
                        help='Number of threads to use (default: all available)')
    opts.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    opts.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}',
                      help='Show program version and exit')

    if len(a) < 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return parser.parse_args(a)


def tigrfam_version_string(accession, max_version=5) -> list[str]:
    """
    Return a list of TIGRFAM version strings for a given accession in descending order from .5 to .1
    TIGR00190 -> ['TIGR00190.5', 'TIGR00190.4', 'TIGR00190.3', 'TIGR00190.2', 'TIGR00190.1']
    """
    accession = accession.split('.')[0] if '.' in accession else accession
    return [
        f'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/{accession}.{i}.HMM' for i in range(max_version, 0, -1)]


def check_file(file: str | Path) -> Path:
    file = Path(file) if isinstance(file, str) else file
    if not file.is_file():
        quit_with_error(f'{file.name} is not a file')
    elif file.stat().st_size == 0:
        quit_with_error(f'{file.name} is empty')
    else:
        return file.absolute()


def is_gzipped(bytes_string: bytes):
    """
    Detects gzipped byte-string
    https://stackoverflow.com/questions/13044562
    """
    return any([bytes_string.startswith(i) for i in [b'\x1f', b'\x8b', b'\x08']])


def load_accessions(accessions: list[str]) -> list['Accession']:
    ok_accessions = []
    for a in accessions:
        if (a := Accession.from_string(a)).urls:
            ok_accessions.append(a)
        else:
            warning(f"Could not determine database for {a}")
    if not ok_accessions:
        quit_with_error("No valid accessions")
    return ok_accessions


def process_accession(acc: 'Accession', outdir: Path = Path.cwd(), extension: str = 'hmm'):
    for url in acc.urls:
        if (hmm := download(url)):
            with open((path := outdir / f"{acc}.{extension}"), 'wb') as f:
                f.write(hmm)
            if path.stat().st_size > 0:  # Check file is not empty
                log(f"Wrote {acc.accession} to {path}")
                break  # Stop trying URLs if one works
            else:
                warning(f"Downloaded {acc.accession} from {url} but file is empty")
                path.unlink()  # Delete empty file and try next URL


def mp_wrapper(arguments):
    """Wrapper for multiprocessing"""
    process_accession(*arguments)


# Classes ---------------------------------------------------------------------
class AccessionError(Exception):
    pass


class Accession:
    def __init__(self, accession: str | None = None, db: str | None = None, urls: list[str] | None = None):
        self.accession = accession or ''
        self.db = db or ''
        self.urls = urls or []

    @classmethod
    def from_string(cls, string: str):
        self = cls(accession=string.strip())
        if self.accession.startswith('PF'):
            self.db = 'PFAM'
            self.urls.append('https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{self.accession}?annotation=hmm')
        elif self.accession.startswith('RF'):
            self.db = 'RFAM'
            self.urls.append(f'https://rfam.org/family/{self.accession}/cm')
        elif self.accession.startswith('cd'):
            self.db = 'CDD'
            self.urls.append(f'https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={self.accession}&seqout')
        elif self.accession.startswith('PTHR'):
            self.db = 'PANTHER'
            self.urls.append(f'http://www.pantherdb.org/panther/exportHmm.jsp?acc={self.accession}')
        elif self.accession.startswith('TIGR') or self.accession.startswith('NF'):
            self.db = 'TIGRFAM'
            self.urls += tigrfam_version_string(self.accession)
        return self

    def __repr__(self):
        return self.accession


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    if len(accessions := load_accessions(args.accessions)) == 1:
        process_accession(args.accessions[0], args.outdir, args.extension)
    else:  # Multiple accessions use multiprocessing
        with Pool(min(len(accessions), args.threads)) as pool:
            pool.imap_unordered(mp_wrapper, [(a, args.outdir, args.extension) for a in accessions])
            pool.close()
            pool.join()

    log("Done!")
