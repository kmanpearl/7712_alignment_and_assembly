# 7712_assembly_and_alignment

This command line application aligns DNA sequence reads to a query sequence, then assembles contigs from all aligned sequences. 
There are two main steps to this program, alignment and assembly. 

1. Alignment: Using a default value of k = 5 the reads and query sequences divided into kmers. 
If an exact match between a read kmer and query sequence kmer is found, the alignment is extended using a dynamic programming algorithm. 
Only alignments above a user inputed threshold (default = ??) are returned. 
1. Assembly: A de Bruijn graph is created from kmers using a user specified k value (default = 10). 
Then a depth-first-search is used to traverse the graph and find all possible paths between start and stop nodes. 

Note that at this time, only the alignment step is completed. 

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

`python main.py -h`

usage: main.py [-h] --q Q --r R [--k K] [--m M] [--mi MI] [--g G] [--t T] [--s S]

Align sequence reads to a query sequence and assemble contigs from aligned reads

optional arguments:
  -h, --help            show this help message and exit
  --k K, -kmer_size K   length of kmers
  --m M, -match_score M
                        alignment score for matching base pairs
  --mi MI, -mismatch_score MI
                        alignment score for non-matching base pairs
  --g G, -gap_score G   alignment score for introducing a gap
  --t T, -score_threshold T
                        minimum normalized score needed to be considered an alignment
  --s S, -save S        if True, save a csv file of the graph

required arguments:
  --q Q, -query_file Q  path to the query FASTA file
  --r R, -read_file R   path to the reads FASTA file

The default length to use when creating kmers is `10`. 
This value can be changed using the optional argument `-kmer_size`. 
The choice of kmer length must be less than the length of the shortest read. 

The default scoring scheme used is +1 for a match and -1 for a mis-match or gap. 
These can be changed using the optional arguments `-match_score`, `-mismatch_score`, and `-gap_score`, respectively. 
This scoring scheme has been chosen for its simplicity; 
however, many other scoring schemes have been developed and may be used instead. 

If `-save` is set to `True` then intermediate output files will be generated. 
These files are:

The default for threshold used is 0.75. 
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

`python scripts/main.py -q "data/QUERY.fasta" -r "data/READS.fatsa" -k 10 -s True`


## Dependencies 