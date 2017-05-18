#!/usr/bin/env python3
'''
  measure gap statistics while building a master gap file
'''

import argparse
import collections
import gzip
import statistics # python 3
import sys

def write_gap(output, gap_info):
    output.write('{}\t{}\t{}\t{:.1f}\n'.format(
            gap_info[0],
            gap_info[1],
            gap_info[1] + len(gap_info[2]),
            statistics.mean(gap_info[2])
        )
    )

def add_stat(stats, gap_info):
    stats['gap_bases'] += len(gap_info[2])
    stats['gap_count'] += 1
    stats['gap_min'] = min(stats['gap_min'], len(gap_info[2]))
    stats['gap_max'] = max(stats['gap_max'], len(gap_info[2]))
    stats['gap_lengths'][len(gap_info[2])] += 1

def find_coverage(coverage_files, output, positions):
    '''
        assume coverage files have the same set of bases
        e.g.
chr1    12099   12227   DDX11L1 1       141
chr1    12099   12227   DDX11L1 2       140
chr1    12099   12227   DDX11L1 3       145
chr1    12099   12227   DDX11L1 4       159
chr1    12099   12227   DDX11L1 5       160
chr1    12099   12227   DDX11L1 6       160
    '''
    # read each coverage file simultaneously
    sys.stderr.write( '{} samples\n'.format(len(coverage_files)))
    fhs = [ gzip.open(coverage, 'rt') for coverage in coverage_files ]
    stats = {'total': 0, 'found': set()}
    sys.stdout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('pos', 'line', 'mean', 'sd', 'n', 'coverages'))
    while True:
        lines = [ fh.readline() for fh in fhs ]
        if lines[0] == '':
            break
        stats['total'] += 1
        fields = [ line.strip('\n').split('\t') for line in lines ]
        # add position
        position = int(fields[0][1]) + int(fields[0][4]) - 1

        if position in positions:
            # calculate mean coverage and sd
            coverages = [ int(field[5]) for field in fields ]
            if len(coverages) > 1:
                mean = statistics.mean(coverages)
                sd = statistics.stdev(coverages)
            else:
                mean = coverages[0]
                sd = 0
            stats['found'].add(position)

            sys.stdout.write('{}\t{}\t{:.2f}\t{:.2f}\t{}\t{}\n'.format(position, stats['total'], mean, sd, len(coverages), ','.join([str(c) for c in coverages])))

        if stats['total'] % 10000 == 0:
            sys.stderr.write( 'processed {} lines. {} gaps.\n'.format(stats['total'], stats['gap_count'] ))

        if len(stats['found']) >= len(positions):
            break

if __name__ == '__main__':
    sys.stderr.write('starting...\n')
    parser = argparse.ArgumentParser(description='Find gaps')
    parser.add_argument('--coverage', nargs='+', required=True, help='coverage files (.gz)')
    parser.add_argument('--position', nargs='+', required=True, help='position to check')
    args = parser.parse_args()
    positions = set([int(p) for p in args.position])
    find_coverage(args.coverage, sys.stdout, positions)
