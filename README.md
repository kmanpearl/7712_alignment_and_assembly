# 7712_assembly_and_alignment

This command line application aligns DNA sequence reads to a query sequence, then assembles contigs from all aligned reads. 
There are two main steps to this program, alignment and assembly. 

1. Alignment: Reads and query sequences are divided into kmers. 
If an exact match between a read kmer and query sequence kmer is found, 
the alignment is extended using a dynamic programming algorithm. 
2. Assembly: A de Bruijn graph is created from k-1mers.
Then a depth-first-search is used to traverse the graph and find all possible paths between start and stop nodes. 

Note that at this time, only the alignment step is completed. 

## Dependencies 

All dependencies are automatically downloaded during the installation process. They are as follows:
```
channels:
  - anaconda
  - conda-forge
dependencies:
  - python=3.9 
  - pip
  - numpy
  - pandas
  - networkx
  - pip:
    - pre-commit
```

## Installation 

For convience, this applications makes use of an anaconda envrionment. 
If you do not already have anaconda installed instructions can be found [here](https://docs.anaconda.com/anaconda/install/). 

First, clone this repository to your current directy be executing the following:

`git clone https://github.com/kmanpearl/7712_alignment_and_assembly.git`


Then create a conda environment from the `envrionment.yml` file. 
This will automatically download all the required dependencies: 

`conda env create -f environment.yml`

To successfully use this application you must first activate the envrionment:

`conda activate align_and_assemble`

## Usage 

To access documentation from the command line run:

`python scripts/main.py -h`

```
usage: main.py [-h] --q Q --r R [--ka KA] [--kb KB] [--m M] [--mi MI] [--g G] [--t T] [--s S]

Align sequence reads to a query sequence and assemble contigs from aligned reads

optional arguments:
  -h, --help            show this help message and exit
  --ka KA, -alignment_kmer KA
                        length of kmers to use for alignment
  --kb KB, -assembly_kmer KB
                        length of kmers to use for assembly
  --m M, -match_score M
                        alignment score for matching base pairs
  --mi MI, -mismatch_score MI
                        alignment score for non-matching base pairs
  --g G, -gap_score G   alignment score for introducing a gap
  --t T, -score_threshold T
                        minimum normalized score needed to be considered an alignment
  --s S, -save S        if True, save intermediate outputs

required arguments:
  --q Q, -query_file Q  path to the query FASTA file
  --r R, -read_file R   path to the reads FASTA file
```

The required arguments are `-query_file`: a FASTA file containing a single query sequence that you wish to align to and `-read_file`: a FASTA file containing all of the sequence reads that you wish to assemble and align against the query. 
The program only accepts DNA sequences reads containing the letters ATGC. 

The default length to use when creating kmers for assembly is `5` for aligning sequence reads. 
This value should be changed to `3` with the argument `-alignment_kmer` if you are running on the files in `test_data`. 
It is recommended to keep this k value small. 
The default value when creating kmers for assembly is `10` but if you are using the `test_data` this value should be changed to `5` using the `-assembly_kmer` argument. 
The choice of kmer length must be less than the length of the shortest read. 

The default scoring scheme used is `+1` for a match and `-1` for a mis-match or gap. 
These can be changed using the optional arguments `-match_score`, `-mismatch_score`, and `-gap_score`, respectively. 
This scoring scheme has been chosen for its simplicity; 
however, many other scoring schemes have been developed and may be used instead. 

The default for threshold used is `0.75`. 
This means that to be considered an alignment a read must have a score (normalized for the length of the read) above 0.75. 
The threshold can be changed using the  `--score_threshold` argument.

The optional argument `-save` is set to `False` by default. 
To save intermediate outputs change this to `True`.
The intermediate output files generated are:


1. `graph.csv`: a two column csv file representing a directed graph.
Column 1 is the source node and column two is the target node. 
1. `fwd_alignment_scores.csv`: normalized alignment scores for each read 
2. `reverse_alignment_scores.csv`: normalized alignment scores for each reversed read 


## Example Use Case 

To run this program on the (very small) test data set with default parameters run the following in the repo directory with the conda envrionment activated:

`python main.py -q "sample_data/query.FASTA.txt" -r "sample_data/reads.FASTA.txt" -ka 3 -kb 5 -s True`

This will run the program with default values and generate all intermediate outputs. 
