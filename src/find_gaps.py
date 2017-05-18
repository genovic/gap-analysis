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

def find_gaps(coverage_files, threshold, sd_offset, output, calculate_stability, max_lines, min_width, gap_filter):
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
    gap_info = None
    stats = {'total': 0, 'gap_bases': 0, 'gap_count': 0, 'gap_min': 1e9, 'gap_max': 0, 'gap_lengths': collections.defaultdict(int)}
    instability_gap = [ 0 ] * len(coverage_files) # counts of stable positions
    instability_nogap = [ 0 ] * len(coverage_files) # counts of stable positions
    while True:
        lines = [ fh.readline() for fh in fhs ]
        if lines[0] == '':
            break
        stats['total'] += 1
        fields = [ line.strip('\n').split('\t') for line in lines ]
        # add position
        position = int(fields[0][1]) + int(fields[0][4]) - 1

        # calculate mean coverage and sd
        coverages = [ int(field[5]) for field in fields ]
        if len(coverages) > 1:
            mean = statistics.mean(coverages)
            sd = statistics.stdev(coverages)
        else:
            mean = coverages[0]
            sd = 0

        if mean + sd * sd_offset < threshold and position not in gap_filter[fields[0][0]]: # it's a gap
            if gap_info is None: # no existing gap - start a new gap
                gap_info = (fields[0][0], position, [mean])
            elif position != gap_info[1] + len(gap_info[2]): # position isn't adjacent - start a new gap
                if len(gap_info[2]) >= min_width: # it's long enough
                    write_gap(output, gap_info)
                    add_stat(stats, gap_info)
                gap_info = (fields[0][0], position, [mean])
            else: # continue existing gap
                gap_info[2].append(mean)

        else: # not a gap
            if gap_info is not None:
                if len(gap_info[2]) >= min_width: # it's long enough
                    write_gap(output, gap_info)
                    add_stat(stats, gap_info)
                gap_info = None

        if calculate_stability:
            was_gap = False
            for idx in range(1, len(coverage_files)):
                mean = statistics.mean(coverages[:idx+1])
                sd = statistics.stdev(coverages[:idx+1])
                is_gap = mean + sd * sd_offset < threshold
                if is_gap != was_gap:
                    if is_gap:
                        instability_gap[idx] += 1
                    else:
                        instability_nogap[idx] += 1
                was_gap = is_gap
        
        if stats['total'] == max_lines:
            break

        if stats['total'] % 10000 == 0:
            sys.stderr.write( 'processed {} lines. {} gaps.\n'.format(stats['total'], stats['gap_count'] ))

    if gap_info is not None:
        write_gap(output, gap_info)
        add_stat(stats, gap_info)

    if calculate_stability:
        stability_output = ' '.join( [ '{:.1f}%'.format(100. * (stab_gap + stab_nogap) / stats['total']) for stab_gap, stab_nogap in zip(instability_gap[1:], instability_nogap[1:]) ] )
        stability_output_gap = ' '.join( [ '{:.1f}%'.format(100. * stab / stats['total']) for stab in instability_gap[1:] ] )
        stability_output_nogap = ' '.join( [ '{:.1f}%'.format(100. * stab / stats['total']) for stab in instability_nogap[1:] ] )
    else:
        stability_output = 'not calculated'
        stability_output_gap = 'not calculated'
        stability_output_nogap = 'not calculated'

    sys.stderr.write('Statistics\n==========\nThreshold: {}\nSD: {}\nBases considered: {}\nBases in gap: {}\nGap count: {}\nMin gap length: {}\nMax gap length: {}\nInstability gap: {}\nInstability no gap: {}\nInstability: {}\nGap lengths: {}\n'.format(
        threshold,
        sd_offset,
        stats['total'], 
        stats['gap_bases'], 
        stats['gap_count'], 
        stats['gap_min'], 
        stats['gap_max'],
        stability_output_gap,
        stability_output_nogap,
        stability_output,
        ' '.join( [ '{}:{}'.format(length, stats['gap_lengths'][length]) for length in sorted(stats['gap_lengths'].keys()) ] )
    ))
    
def build_filter(in_fh):
    sys.stderr.write('building master gap filter...')
    result = collections.defaultdict(set)
    for idx, line in enumerate(in_fh):
        fields = line.strip('\n').split('\t')
        if len(fields) >= 3:
            for pos in range(int(fields[1]), int(fields[2])):
                result[fields[0]].add(pos)
        if idx % 10000 == 0:
            sys.stderr.write( 'processed {} lines.\n'.format(idx))
    sys.stderr.write('building master gap filter: done')
    return result

if __name__ == '__main__':
    sys.stderr.write('starting...\n')
    parser = argparse.ArgumentParser(description='Find gaps')
    parser.add_argument('--coverage', nargs='+', required=True, help='coverage files (.gz)')
    parser.add_argument('--threshold', type=int, default=15, help='consider a gap if mean coverage falls below this')
    parser.add_argument('--sd', type=float, default=-1, help='offset mean by this proportion of the standard deviation when considering if a base is a gap')
    parser.add_argument('--min_width', type=float, default=1, help='only report gaps if the length is at least this')
    #parser.add_argument('--median', type=bool, action='store_true', default=False, help='use median statistics instead of mean')
    parser.add_argument('--stability', action='store_true', default=False, help='calculate stability as more samples added')
    parser.add_argument('--max_lines', type=int, default=1e9, help='stop after reading this many lines')
    parser.add_argument('--filter', required=False, help='bed file specifying parts of the genome to ignore')
    args = parser.parse_args()
    if args.filter:
        gap_filter = build_filter(open(args.filter, 'r'))
    else:
        gap_filter = collections.defaultdict(set)
    find_gaps(args.coverage, args.threshold, args.sd, sys.stdout, args.stability, args.max_lines, args.min_width, gap_filter)
