# 7712_assembly_and_alignment

This command line program assembles DNA sequence reads from a FASTA file. 
The assembled contigs are then aligned againts a query sequence.
There are two main components to this program.

1. **Assembly:** A de Bruijn graph is created with k-1mers as nodes and kmers as directed edges.
Then a depth-first-search is used to traverse the graph and find all possible paths between start and stop nodes.
Each path is used to assemble contigs 
1. **Alignment:** Reads and assembled contigs are divided into kmers. 
If an exact match between a read kmer and query sequence kmer is found, 
the alignment is extended using a dynamic programming algorithm. 
All scores above a user specified threshold are considered an alignment. 




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
  - pip:
    - pre-commit
```

## Installation 

For convience, this applications makes use of an anaconda envrionment. 
If you do not already have anaconda installed instructions can be found [here](https://docs.anaconda.com/anaconda/install/). 

First, clone this repository to your current directy be executing the following from the terminal:

`git clone https://github.com/kmanpearl/7712_alignment_and_assembly.git`


Then create a conda environment from the `envrionment.yml` file. 
This will automatically download all the required dependencies: 

`conda env create -f environment.yml`

To successfully use this application you must first activate the envrionment:

`conda activate align_and_assemble`


Currently the file pathing is inflexible.
The output directory named `output` must exist.
If this directory is preexisting, any files in it with the same name as ouput files will be overwritten.
To create the output directory:

`mkdir output`

## Usage 

To access documentation from the command line run:

`python scripts/main.py -h`

```
usage: main.py [-h] --q Q --r R [--k K] [--m M] [--mi MI] [--g G] [--t T] [--s S]

Assemble sequence reads and align to a query

optional arguments:
  -h, --help            show this help message and exit
  --k K, -kmer_size K   length of kmers: must be shorter than the shortest read
  --m M, -match_score M
                        alignment score for matching base pairs
  --mi MI, -mismatch_score MI
                        alignment score for non-matching base pairs
  --g G, -gap_score G   alignment score for introducing a gap
  --t T, -score_threshold T
                        minimum normalized score needed to be considered an alignment:
                        value must be between 0-1
  --s S, -save S        if True, save intermediate outputs

required arguments:
  --q Q, -query_file Q  path to the query FASTA file
  --r R, -read_file R   path to the reads FASTA file
```

The required arguments are `-query_file`: a FASTA file containing a single query sequence that you wish to align to and `-read_file`: a FASTA file containing all of the sequence reads that you wish to assemble and align against the query. 
The reads file must contain more than one sequence.
The program only accepts DNA sequences reads containing the letters ATGC.
Any other characters, such as those used to represent ambiguos base pairs, will raise an exception.  


The default length to use when creating kmers is `15`. 
This value should be changed to `5` with the argument `-kmer_size` if you are running on the files in `test_data`. 
It is recommended to keep this k value small. 
If the length of kmer is shorter than the smallest sequence read it will raise an exception. 

The default scoring scheme used is `+1` for a match and `-1` for a mis-match or gap. 
These can be changed using the optional arguments `-match_score`, `-mismatch_score`, and `-gap_score`, respectively. 
This scoring scheme has been chosen for its simplicity; 
however, many other scoring schemes have been developed and may be used instead. 

Alignment scores are created by dividing the best score in the alignment matrix by the length of the read. 
Scores range from 0-1 where 0 represents no matching positions and 1 represents a perfect match. 
The default for the score threshold used is `0.5`. 
This means that to be considered an alignment, a read must have a score above 0.5. 
The threshold can be changed using the  `-score_threshold` argument.

The optional argument `-save` is set to `False` by default. 
To save intermediate outputs change this to `True`.
The intermediate output files generated are `adjacency_matrix.csv` and `alignment_scores.csv`.



## Output

All outputs will be generated in the directory `output` which must exist.
If the directory does not exist, please make one by running `mkdir output` from the terminal when in the project directory. 

The output files generated are:

1. `ALLELES.fasta`: a fasta file containing the sequence of the longest contig that aligned to the query
2. `ALLELES.aln`: a csv containing information about the sequence alignment. 
The columns are as follows:
`sseqid` the sequence id of the aligned read
`cseqid` the contig id 
`sstart` the starting position where the read aligned to the contig
`ssend` the ending position where the read aligned to the contig
`qstart` the starting position in the contig that the read aligned too
`qend` the ending position in the contig that the read aligned too

The optional output files generated are: 

1. `adjacency_matrix.csv`: a table containing the adjacency matrix made from the graph. 
Rows are source nodes and columns target nodes. 
1. `alignment_scores.csv`: normalized alignment scores for each read (only reported for reads above the user specified threshold)

## Example Use Case 

To run this program on the (very small) test data set run the following in the cloned repo directory, with the conda envrionment activated:

`python main.py -q "sample_data/query.FASTA.txt" -r "sample_data/reads.FASTA.txt" -k 5 -s True`

