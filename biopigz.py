#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
biopigz
-----------------
Requires:
- Python >=3.9
-----------------
Copyright 2022 Tom Stanton (tomdstanton@gmail.com)
"""
# ............... Base imports ............... #
from argparse import ArgumentParser, RawTextHelpFormatter
import logging
from sys import stdin, stdout, exit, argv
import zlib
from multiprocessing.dummy import Pool
# import cython
from os import cpu_count
from pathlib import Path
from time import time
from uuid import uuid4

LOGGER = logging.getLogger(__name__)
__version__ = '0.0.1'


class InputFile:
    def __init__(self, path, _args):
        self.path = Path(path) if path != '-' else '-'
        self.out_file = None
        self.stream = self.read()
        self.converted_stream = None
        self.stdout = _args.stdout
        self.keep = _args.keep
        self.level = _args.level

    def read(self):
        if self.path == '-':
            return stdin.buffer.read()
        else:
            with open(self.path, 'rb') as fh:
                return fh.read()

    def write(self):
        if self.stdout:
            stdout.write(self.converted_stream.decode())
        else:
            start = time()
            with open(self.out_file, 'wb') as fh:
                fh.write(self.converted_stream)
            LOGGER.info(f'Wrote {self.out_file} in {time() - start:.2f} seconds')
            if not self.keep and self.path != '-':
                self.path.unlink()
            LOGGER.info(f'Deleted {self.path}')


class CompressedFile(InputFile):
    def __init__(self, path, _args):
        super().__init__(path, _args)
        self.decoder = zlib.decompressobj(wbits=zlib.MAX_WBITS | 32)
        self.out_file = path.with_suffix('') if path != '-' else Path().cwd() / f'{uuid4()}.file'

    def convert(self):
        start = time()
        self.converted_stream = self.decoder.decompress(self.stream) + self.decoder.flush()
        LOGGER.info(f'{self.path} decompressed '
                    f'{100 * len(self.converted_stream) / len(self.stream):.2f}% in {time() - start:.2f} seconds')


class DecompressedFile(InputFile):
    def __init__(self, path, _args):
        super().__init__(path, _args)
        self.decoder = zlib.compressobj(level=args.level)
        self.out_file = self.path.with_suffix(self.path.suffix + '.gz') if path != '-' \
            else Path().cwd() / f'{uuid4()}.file.gz'

    def convert(self):
        start = time()
        self.converted_stream = self.decoder.compress(self.stream) + self.decoder.flush()
        LOGGER.info(f'{self.path} compressed '
                    f'{100 * len(self.converted_stream) / len(self.stream):.2f}% in {time() - start:.2f} seconds')


def parse_args(_args):
    parser = ArgumentParser(
        add_help=False, prog='biopigz',
        description='A simple Python gzip tool, designed for processing multiple medium-sized files in parallel',
        formatter_class=RawTextHelpFormatter,
        usage='biopigz [options] [files ...]',
        epilog=f" - Processes files in place, adding/removing the '.gz' suffix.\n"
               f" - If stdin and not stdout, processed output will be written to [uuid.file(.gz)]\n"
               f" - Multiple files will be processed simultaneously using --threads\n"
               #f" - If stdin, chunks will be processed simultaneously using --threads \n"
               )

    input_opts = parser.add_argument_group('Input')
    input_opts.add_argument('files', nargs='*', help='Files to process, or stdin if none provided',
                        type=lambda p: Path(p).absolute())

    compression_opts = parser.add_argument_group('Compression options')
    compression_opts.add_argument('-d', '--decompress', default=False, action='store_true',
                                  help="Decompress the compressed input (files must have the '.gz' suffix)")
    compression_opts.add_argument('-c', '--stdout', default=False, action='store_true',
                                  help="Write all processed output to stdout (won't delete input files)")
    compression_opts.add_argument('-k', '--keep', default=False, action='store_true',
                                  help="Do not delete original file after processing")
    compression_opts.add_argument('-l', '--level', type=int, metavar='0-9', choices=range(0, 10), default=6,
                                  help="Compression level")

    other_opts = parser.add_argument_group('Other options')
    other_opts.add_argument('-t', '--threads', default=(x := cpu_count()), metavar='N',
                            type=lambda t: min(x, int(t)), help=f'Number of threads, default: {x}')
    other_opts.add_argument('--version', action='version', version=f'{__name__} v{__version__}',
                            help="Show version number and exit")
    other_opts.add_argument('-v', '--verbosity', action="count", default=0, help="Verbosity, increases with -vv")
    other_opts.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    return check_args(parser.parse_args(_args))


def check_args(_args):
    if _args.verbosity == 1:
        _args.verbosity = logging.INFO
    elif _args.verbosity > 1:
        _args.verbosity = logging.DEBUG
    else:
        _args.verbosity = logging.WARNING
    logging.basicConfig(format="[%(asctime)s %(levelname)s %(funcName)s] %(message)s", datefmt='%H:%M:%S',
                        level=_args.verbosity)
    if _args.files:
        _args.files = [p for p in _args.files if check_file(p, _args.decompress)]
        if not _args.files:
            quit_with_error('No valid files found')
    LOGGER.info(' '.join(f'{k}={v}' for k, v in vars(_args).items() if k != 'files'))
    return _args


def check_extension(path: Path, d=False):
    return path.suffix == '.gz' if d else path.suffix != '.gz'


def check_file(_file: Path, d=False):
    if not _file.is_file():
        LOGGER.warning(format_str(f'{_file.name} does not exist, skipping', '93'))
        return None
    elif _file.stat().st_size < 1:
        LOGGER.warning(format_str(f'{_file.name} is empty, skipping', '93'))
        return False
    elif not check_extension(_file, d):
        LOGGER.warning(format_str(f'{_file.name} must not end in ".gz", or end in ".gz" for -d, skipping', '93'))
        return False
    else:
        return True


def process_file(file, _args):
    if file != '-' and not check_file(file, _args.decompress):
        return None
    else:
        infile = CompressedFile(file, _args) if args.decompress else DecompressedFile(file, _args)
        infile.convert()
        infile.write()


def format_str(string, format_number):
    return f"\033[{format_number}m{string}\033[0m"


def quit_with_error(message):
    """Displays the given message and ends the program's execution."""
    LOGGER.error(format_str(message, '91'))
    exit(1)


def task_wrapper(arguments):
    return process_file(*arguments)


if __name__ == '__main__':
    args = parse_args(argv[1:])
    # args = parse_args(['ERR4920592_1.fastq.gz', '-d'])

    if not args.files:
        if stdin.isatty():
            quit_with_error('No input provided, exiting')
        else:
            process_file('-', args)

    else:
        with Pool(min(len(args.files), args.threads)) as pool:
            pool.imap_unordered(task_wrapper, [(file, args) for file in args.files])
            pool.close()
            pool.join()

    exit(0)
