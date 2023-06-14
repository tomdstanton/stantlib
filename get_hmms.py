#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'Tom Stanton'
__title__ = 'Get HMMs'
__author_email__ = 'tomdstanton@gmail.com'
__description__ = 'Download HMMs from public databases with accessions'
__license__ = 'gpl-3.0'
__version__ = '0.0.1'

# ................ Python Imports ................ #
from pathlib import Path
import requests
from requests.exceptions import HTTPError
import sys
from argparse import RawTextHelpFormatter, ArgumentParser
import logging
from io import BytesIO
import multiprocessing.pool
from time import sleep
# ................ Dependency Imports ................ #
import pyhmmer
# ................ Globals ................ #
LOGGER = logging.getLogger(__name__)  # https://stackoverflow.com/a/45065655/10771025


class Accession:
    def __init__(self, accession: str, out_dir: Path):
        self.accession = accession.strip()
        self.db, self.url = determine_db(self.accession)
        self.out_path = Path(out_dir / self.accession.replace(":", "-"))
        self.hmm_path = Path(str(self.out_path) + ".hmm")
        self.hmm = None

    def __repr__(self):
        return self.accession

    def get(self):
        if self.url:
            download_string = self.download_hmm()
            if download_string:
                self.hmm = self.check_hmm(download_string)
                if not self.hmm:
                    LOGGER.warning(format_str(f"Could not generate a proper HMM for {self.accession}", '93'))
                else:
                    self.write_hmm()

    def check_hmm(self, download_string: str):
        LOGGER.info(f"Checking {self.accession}")
        try:
            with BytesIO() as fh:
                fh.write(download_string.encode())
                fh.seek(0)
                with pyhmmer.plan7.HMMFile(fh) as hmm_file:
                    hmm = hmm_file.read()

        except Exception as e:
            LOGGER.warning(format_str(str(e), '93'))
            hmm = build_hmm_from_msa(download_string, self.accession)
        return hmm

    def write_hmm(self):
        try:
            with open(self.hmm_path, 'wb') as fh:
                self.hmm.write(fh)
        except Exception as e:
            LOGGER.warning(format_str(str(e), '93'))

        if self.hmm_path.stat().st_size > 0:
            LOGGER.info(f"HMM successfully written to {self.hmm_path}")

    def download_hmm(self, retries=3):
        LOGGER.info(f"Downloading {self.accession} from {self.db}")
        hmm_string = b''
        for n in range(retries):
            try:
                with requests.get(self.url, stream=True) as request:
                    request.raise_for_status()
                    for chunk in request.iter_content(chunk_size=2 ** 20):
                        hmm_string += chunk
                break
            except HTTPError as exc:
                code = exc.response.status_code
                if code in [429, 500, 502, 503, 504]:
                    sleep(n)  # retry after n seconds
                    continue
                raise
        if not hmm_string:
            LOGGER.warning(format_str(f"Failed to download HMM from {self.url}, check accession", '93'))
            return None

        LOGGER.debug(f"Downloaded {len(hmm_string)} bytes")
        return hmm_string.decode()


def parse_args(arguments):
    parser = ArgumentParser(formatter_class=RawTextHelpFormatter,
                            description=__description__)
    parser.add_argument('search_terms', nargs='+', type=str)
    parser.add_argument('-o', '--outdir', type=lambda p: check_dir(parser, Path(p).resolve()),
                        default=Path.cwd(),
                        help="Directory to save HMMs")
    return parser.parse_args(arguments)


def check_dir(parser: 'argparse.ArgumentParser', dirpath) -> Path:
    try:
        dirpath.mkdir(parents=True, exist_ok=True)
        return dirpath
    except Exception as e:
        parser.error(str(e))


def build_hmm_from_msa(fasta_string: str, name: str):
    alphabet = pyhmmer.easel.Alphabet.amino()
    seqs = []

    fasta_string = "\n".join(
        [i for i in fasta_string.split('\n') if i[0].upper() in alphabet.symbols + '>']
    )
    try:
        for record in fasta_string.split('\n>'):
            seqs += [
                pyhmmer.easel.TextSequence(
                    name=record.split('\n')[0].split(' ')[0].replace('>', '').encode(),
                    sequence=''.join(record.split('\n')[1:])
                )
            ]
    except Exception as e:
        LOGGER.warning(format_str(str(e), '93'))
        return None

    try:
        msa = pyhmmer.easel.TextMSA(
            name=name.encode(),
            sequences=seqs
        ).digitize(alphabet)
        LOGGER.info(f"{name} successfully parsed as alignment")
    except Exception as e:
        LOGGER.warning(format_str(str(e), '93'))
        return None

    try:
        builder = pyhmmer.plan7.Builder(alphabet)
        background = pyhmmer.plan7.Background(alphabet)
        hmm, _, _ = builder.build_msa(msa, background)
    except Exception as e:
        LOGGER.warning(format_str(str(e), '93'))
        return None

    return hmm


def determine_db(accession: str):
    if accession.startswith('PF'):
        return 'PFAM', f'https://pfam-legacy.xfam.org/family/{accession}/hmm'
    # elif accession.startswith('RF'):
    #     return 'RFAM', f'https://rfam.org/family/{accession}/cm'
    elif accession.startswith('cd'):
        return 'CDD', f'https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={accession}&seqout'
    elif accession.startswith('PTHR'):
        return 'PANTHER', f'http://www.pantherdb.org/panther/exportHmm.jsp?acc={accession}'
    elif accession.startswith('TIGR'):
        return 'TIGRFAM', f'https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/{accession}.HMM'
    else:
        LOGGER.warning(format_str(f"Could not determine database for {accession}", '93'))
        return None, None


def format_str(string: str, format_number: str) -> str:
    """Add colour to string using a format string number"""
    return f"\033[{format_number}m{string}\033[0m"


def quit_with_error(message: str):
    """Displays the given message and ends the program's execution."""
    LOGGER.error(format_str(f"ERROR: {message}", '91'))
    sys.exit(1)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    line_width = 70  # Nice general terminal width
    logging.basicConfig(format="[%(asctime)s %(levelname)s %(funcName)s] %(message)s",
                        datefmt='%H:%M:%S', level=logging.INFO)
    LOGGER.info(f' {Path(__file__).stem} {__version__} '.center(line_width, '='))
    LOGGER.info(f'Platform: {sys.platform}')
    LOGGER.info(f'Python: {".".join([str(i) for i in sys.version_info[:3]])}')
    LOGGER.info(f'Current working directory: {Path.cwd()}')

    def f(acc: Accession):
        return acc.get()

    with multiprocessing.pool.ThreadPool(multiprocessing.cpu_count()) as pool:
        results = pool.map(f, [Accession(i, args.outdir) for i in args.search_terms])
