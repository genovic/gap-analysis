#!/usr/bin/env python
'''
  filters a bed file with a list of genes
  e.g.
  python filter_genes.py gene_list.txt < bed_file > filtered_bed_file
'''

import sys


def main():
    genes = set()
    for line in open(sys.argv[1], 'r'):
        genes.add(line.strip('\n').split('\t')[0])

    sys.stderr.write('{} genes\n'.format(len(genes)))

    included = 0
    total = 0
    for total, line in enumerate(sys.stdin):
        fields = line.strip('\n').split('\t')
        if fields[3] in genes:
            sys.stdout.write(line)
            included += 1

        if total % 1000000 == 0:
            sys.stderr.write('included {} of {}\n'.format(included, total))

    sys.stderr.write('included {} of {}\n'.format(included, total))


if __name__ == '__main__':
    main()
