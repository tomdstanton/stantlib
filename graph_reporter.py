#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
graph_reporter
-----------------
Uses Bandage info to generate a formatted report for
a number of assembly graphs.
-----------------
Requires:
- Python >=3.9
- Bandage
-----------------
Tom Stanton, 2021
"""

from pathlib import Path
from argparse import ArgumentParser
import logging
import os
import subprocess
import sys

__version__ = '0.0.1'
__author__ = "Tom Stanton"
__maintainer__ = "Tom Stanton"
__email__ = "tomdstanton@gmail.com"
__status__ = "Development"


def main():
    # ............... Start Program ...............
    args = parse_args()
    logger = logging.getLogger('root')
    logging.basicConfig(format="%(asctime)s | %(levelname)s | %(funcName)s: %(message)s", datefmt="%H:%M:%S")
    logger.setLevel('INFO')

    out_tab = ''
    header = 'Sample\tNode count\tEdge count\tSmallest edge overlap (bp)\tLargest edge overlap (bp)\t' \
             'Total length (bp)\tTotal length no overlaps (bp)\tDead ends\tPercentage dead ends\t' \
             'Connected components\tLargest component (bp)\tTotal length orphaned nodes (bp)\tN50 (bp)\t' \
             'Shortest node (bp)\tLower quartile node (bp)\tMedian node (bp)\tUpper quartile node (bp)\t' \
             'Longest node (bp)\tMedian depth\tEstimated sequence length (bp)'

    num_inputs = len(args.input)
    logger.info(f'Generating report for {num_inputs} assembly graphs')
    counter = 1
    for gfa_file in args.input:
        print(f'[{counter}/{num_inputs}]', end="\r", flush=True, file=sys.stderr)
        counter += 1

        if not os.path.isfile(gfa_file):
            logger.error(f'{gfa_file} is not a valid file')
            continue
        filename = os.path.splitext(os.path.basename(str(gfa_file)))[0]
        info = parse_Bandage_info(gfa_file)
        if info:
            out_tab += f'{filename}\t{info}'

    if out_tab != '':
        if args.out is not None:
            outfile = args.out
            if outfile == '':
                outfile = os.path.join(os.getcwd(), 'graph_report.tab')
            with open(outfile, 'wt') as out:
                out.write(header + '\n' + out_tab.strip())
            logger.info(f'Written report to: {outfile}')
        else:
            print(header + '\n' + out_tab.strip())
    else:
        logger.info('Nothing to output')


def parse_Bandage_info(gfa_path):
    logger = logging.getLogger('root')
    bandage_out = run_Bandage_info(gfa_path)
    formatted_string = ''
    if bandage_out:
        for line in bandage_out.split('\n'):
            formatted_string += line.split(':')[1].strip() + '\t'
        return formatted_string.strip() + '\n'  # removes trailing tab, add newline
    else:
        logger.error(f'No result for: {gfa_path}')


def run_Bandage_info(gfa_path):
    logger = logging.getLogger('root')
    cmd = ['Bandage', 'info', gfa_path]
    child = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = child.communicate()
    if out:
        return out.decode('utf-8').strip()
    if err:
        logger.warning(err.decode("utf-8").strip())


def parse_args():
    parser = ArgumentParser(add_help=False, usage="graph_reporter.py *.gfa")
    parser.print_usage = parser.print_help
    parser.add_argument('input', nargs='+', type=Path)
    parser.add_argument('-o', '--out', nargs='?', type=str, const='',
                        help='Write to file, if no path specified write graph_report.tab in CWD')
    return parser.parse_args()


if __name__ == '__main__':
    main()
