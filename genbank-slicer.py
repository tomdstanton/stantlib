#!/usr/bin/env python3
# -*- coding: utf-8 -*-



from Bio import SeqIO
import argparse
import sys
import os
import re



def main():
    args = parse_args(sys.argv[1:])  # Uncomment below for testing
    #args = parse_args('/Users/tom/Bioinformatics/python_scripts/gff-slicer/GCF_020526025.1.gbk -if TonB'.split())

    start, end = None, None
    if args.slice:
        try:
            start, end = args.slice.split(':')
            start = int(start) if start.isdigit() else str(start)
            end = int(end) if end.isdigit() else str(end)
            if start == end:
                quit_with_error(f'Start slice point {start} is the same as end slice point {end}, '
                                f'please provide two different slice points')
        except Exception as e:
            quit_with_error(f'{e}\n{args.slice} not a valid slice format')


    for gbk in args.genbank:
        # \\\ Check input \\\ #
        if gbk == '-':
            gbk = sys.stdin
        else:
            if not os.path.isfile(gbk):
                print(f'[FILE]: {gbk} is not valid', file=sys.stderr)
                continue

        # \\\ Iterate over records \\\ #
        for rec in SeqIO.parse(gbk, 'gff'):
            if args.include_records:
                x = [i for i in args.include_records if i in [rec.id, rec.name]]
                if not x:
                    print(f"[EXCLUDE]: {rec.id} ({rec.description})", file=sys.stderr)
                    continue
            if args.exclude_records:
                x = [i for i in args.exclude_records if i in [rec.id, rec.name]]
                if x:
                    print(f"[EXCLUDE]: {rec.id} ({rec.description})", file=sys.stderr)
                    args.exclude_records = [i for i in x if i not in args.exclude_records]
                    continue

            print(f"[RECORD]: {rec.id} ({rec.description})", file=sys.stderr)

            # \\\ Init the slicing points \\\ #
            slice_points = []
            # \\\ If args are integers, set the slicing points \\\ #
            if type(start) == int:  # Just use int arguments as slice points
                slice_points.append(start)
            if type(end) == int:  # Just use int arguments as slice points
                slice_points.append(end)
            # \\\ If not, search feature qualifiers for start and end markers \\\ #
            matched_start_marker, matched_end_marker = 0, 0

            if type(start) == str or type(end) == str:  # if either terms are strings, iterate through
                for feature in rec.features:
                    quals = [feature.qualifiers[x] for x in [i for i in args.qualifiers if i in feature.qualifiers.keys()]]
                    quals = set([item for sublist in quals for item in sublist])  # Flatten the qualifiers list of single-element lists

                    # if args.include_features:
                    #     x = [i for i in args.include_features if i in quals]
                    #     if not x:
                    #         print(f"[EXCLUDE]: feature at "
                    #               f"{int(feature.location.start)}:{int(feature.location.end)} ", file=sys.stderr)
                    #         rec.features.pop(feature)
                    #
                    # if args.exclude_features:
                    #     x = [i for i in args.exclude_features if i in quals]
                    #     if x:
                    #         print(f"[EXCLUDE]: feature at "
                    #               f"{int(feature.location.start)}:{int(feature.location.end)} ", file=sys.stderr)
                    #         rec.features.pop(feature)

                    if type(start) == str:
                        start_pattern = re.compile(start, re.IGNORECASE)
                        if any((match := start_pattern.match(x)) for x in quals):
                            if int(feature.location.start) and int(feature.location.end) not in slice_points:
                                slice_points += int(feature.location.start), int(feature.location.end)
                                print(f"[MATCH]: Start marker {start} with {match.group(0)} "
                                      f"({int(feature.location.start)}:{int(feature.location.end)})", file=sys.stderr)
                                matched_start_marker += 1

                    if type(end) == str:
                        end_pattern = re.compile(end, re.IGNORECASE)
                        if any((match := end_pattern.match(x)) for x in quals):
                            if int(feature.location.start) and int(feature.location.end) not in slice_points:
                                slice_points += int(feature.location.start), int(feature.location.end)
                                print(f"[MATCH]: End marker {end} with {match.group(0)} "
                                      f"({int(feature.location.start)}:{int(feature.location.end)})", file=sys.stderr)
                                matched_end_marker += 1

            # \\\ Perform final checks before slicing the record \\\ #
            if args.slice:
                if type(start) == str and matched_start_marker == 0:
                    print(f"[SKIP RECORD]: Could not find a slice point for {start}, try a different marker", file=sys.stderr)
                    continue

                if type(end) == str and matched_end_marker == 0:
                    print(f"[SKIP RECORD]: Could not find a slice point for {end}, try a different marker", file=sys.stderr)
                    continue

                if matched_start_marker > 1:
                    print(f"[SKIP RECORD]: More than one potential slice point for {start}, use a more specific marker", file=sys.stderr)
                    continue

                if matched_end_marker > 1:
                    print(f"[SKIP RECORD]: More than one potential slice point for {end}, use a more specific marker", file=sys.stderr)
                    continue

                start_slice, end_slice = min(slice_points), max(slice_points)

                if end_slice > len(rec):
                    print(f"[SKIP RECORD]: {end} ({end_slice}) is > the length of the record ({len(rec)})", file=sys.stderr)
                    continue

                if start_slice == end_slice:
                    print(f"[SKIP RECORD]: Start slice point {start} is == end slice point {end}", file=sys.stderr)
                    continue

                # \\\ Perform final slice and write to STDOUT \\\ #
                print(f"[SLICING]: Start slice point {start_slice}, end slice point {end_slice}", file=sys.stderr)
                SeqIO.write(rec[start_slice:end_slice], sys.stdout, args.outfmt)
            else:
                SeqIO.write(rec, sys.stdout, args.outfmt)


def parse_args(args):
    parser = argparse.ArgumentParser(description='gff-slicer',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gff', nargs='+', help='Genbank file(s) or stream', default='-', metavar='<STDIN>')
    parser.add_argument('-s', '--slice', nargs='?',
                        help=f'Slice range, formatted as start:end.\n'
                             f'This can be a mix of integers and strings.\n'
                             f'If an integer is included, the slicing will be performed at that exact location.\n'
                             f'If a string is included (e.g. a gene name), various feature attributes will be searched '
                             f'for that particular string (case insensitive).\n'
                             f'An error will be thrown if there are multiple slice points found for a particular string.')
    # parser.add_argument('-if', '--include_features', nargs='*', help='Space-separated list of features to include')
    # parser.add_argument('-ef', '--exclude_features', nargs='*', help='Space-separated list of features to exclude')
    parser.add_argument('-ir', '--include_records', nargs='*', help='Space-separated list of records to include')
    parser.add_argument('-er', '--exclude_records', nargs='*', help='Space-separated list of records to exclude')
    parser.add_argument('-o', '--outfmt', choices=['gff', 'fasta'], help='<STOUT format>', default='gff')
    parser.add_argument('-q', '--qualifiers', nargs='*', help='Genbank qualifiers to search for string slice points',
                        default=['gene', 'locus_tag', 'old_locus_tag', 'note', 'product', 'protein_id'])
    parser.add_argument('--version', action='version', version=f'{__title__} {__version__}',
                        help="Show version number and exit")
    return parser.parse_args(args)


if __name__ == '__main__':
    main()
