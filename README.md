# Gap assessment

## Installation

* Run:
    ```
    git clone git@github.com:genovic/gap-analysis.git
    cd gap-analysis
    python3 -m venv venv
    source venv/bin/activate
    pip install .
    ```
    
* The three commands will now be available:
   * `find_gaps`
   * `find_coverage`
   * `filter_genes`

## Usage

### Find Gaps

```
usage: find_gaps [-h] --coverage COVERAGE [COVERAGE ...]
                 [--threshold THRESHOLD] [--sd SD] [--min_width MIN_WIDTH]
                 [--stability] [--max_lines MAX_LINES] [--filter FILTER]

Find gaps

optional arguments:
  -h, --help            show this help message and exit
  --coverage COVERAGE [COVERAGE ...]
                        coverage files (.gz)
  --threshold THRESHOLD
                        consider a gap if mean coverage falls below this
  --sd SD               offset mean by this proportion of the standard
                        deviation when considering if a base is a gap
  --min_width MIN_WIDTH
                        only report gaps if the length is at least this
  --stability           calculate stability as more samples added
  --max_lines MAX_LINES
                        stop after reading this many lines
  --filter FILTER       bed file specifying parts of the genome to ignore
```

## Find Coverage

```
usage: find_coverage [-h] --coverage COVERAGE [COVERAGE ...] --position
                     POSITION [POSITION ...]

Find gaps

optional arguments:
  -h, --help            show this help message and exit
  --coverage COVERAGE [COVERAGE ...]
                        coverage files (.gz)
  --position POSITION [POSITION ...]
                        position to check
```

## Filter Genes
```
usage: filter_genes gene_list.txt < bed_file > filtered_bed_file
```
