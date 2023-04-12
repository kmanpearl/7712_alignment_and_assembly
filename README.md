# 7712_assembly_and_alignment

This command line program assembles DNA sequence reads from a FASTA file. 
The assembled contigs are then aligned againts a query sequence.
There are two main components to this program are:

1. **Assembly:** A de Bruijn graph is created with k-1mers as nodes and kmers as directed edges.
Then a depth-first-search is used to traverse the graph and find all possible paths between start and stop nodes.
Each path is used to assemble contigs.
2. **Alignment:** The query sequence and assembled contigs are divided into kmers. 
If an exact match between a read kmer and query sequence kmer is found, 
the alignment is extended using a dynamic programming algorithm. 
All alignment scores above a user specified threshold are considered a true alignment. 




## Dependencies 

All dependencies are automatically downloaded during the installation process. They are as follows:

```
  - python=3.9 
  - pip
  - numpy
  - pandas
  - pre-commit
```

## Installation 

For convience and reproducability, this applications makes use of an anaconda envrionment. 
If you do not already have anaconda installed instructions can be found [here](https://docs.anaconda.com/anaconda/install/). 

First, clone this repository to your current directy be executing the following from the terminal:

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
usage: main.py [-h] --q Q --r R --o O [--k K] [--m M] [--mi MI] [--g G]
               [--t T] [--s S]

Assemble sequence reads and align to a query

required arguments:
  --q Q, -query_file Q  path to the query FASTA file
  --r R, -read_file R   path to the reads FASTA file
  --o O, -output_directory O
                        directory to store all generated output files

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

```

The required arguments are 

`-query_file`: a FASTA file containing a single query sequence that you wish to align to and 

`-read_file`: a FASTA file containing all of the sequence reads that you wish to assemble and align against the query. 
`-output_directory`: directory name to store all outputs generated. 
The directory must exist already. 
To create the directory run `mkdir <output_directory>` in the terminal window.


The query file may only contain one sequence, while the reads file must contain more than one sequence.
The program only accepts DNA sequences reads containing the letters ATGC.
Any other characters, such as those used to represent ambiguos base pairs, will raise an exception.  

The default length to use when creating kmers is `30`. 
This value can be changed with the argument `-kmer_size`.

If the length of kmer is shorter than the smallest sequence read it will raise an exception. 
Choice of kmer length should reflect the similarity between sequences.
For sequences with more variation, a smaller k value will be more likely to capture all alignments.
For sequences with little variation, a larger k value will be more stringent in its selection criteria. 
This reduces the number of sequences selected for alignment, thus the runtime of the program.

The default scoring scheme used is `+1` for a match and `-1` for a mis-match or gap. 
These can be changed using the optional arguments `-match_score`, `-mismatch_score`, and `-gap_score`, respectively. 
This scoring scheme has been chosen for its simplicity; 
however, many other scoring schemes have been developed and may be used instead. 

Alignment scores are created by dividing the best score in the alignment matrix by the length of the read. 
Scores range from 0-1 where 0 represents no matching positions and 1 represents a perfect match. 
The default for the score threshold used is `0.75`. 
This means that to be considered an alignment, a read must be more than 75% similar to the query sequence.
The threshold can be changed using the  `-score_threshold` argument.
Like kmer length, the score threshold shouold reflect the similarity between sequences.
The more related the sequences are the higher the threshold that should be used.

The optional argument `-save` is set to `False` by default. 
To save intermediate outputs change this to `True`.
The intermediate output files generated are `adjacency_matrix.csv` and `alignment_scores.csv`.



## Output

All outputs will be generated in the directory specified by `-output_dir`.
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

To run this program on the (very small) test data set and generate intermediate outputs, run the following in the cloned repo: 

`mkdir output`

`conda activate align_and_assemble`

`python main.py -q "sample_data/fake_QUERY.fasta" -r "sample_data/fake_READS.fasta" -o output -s True`

