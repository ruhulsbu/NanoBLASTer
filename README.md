## NanoBLASTer: Basic Local Alignment and Search Tool for Oxford Nanopore Long Sequences  
**__Current Version: 0.16__**  
Release Date: October 13, 2015  
Platform: Linux x64 system

**Nanopore sequencing data**: [Oxford Nanopore Data for Reference and Reads](http://schatzlab.cshl.edu/data/nanocorr/)

**Please contact moamin@cs.stonybrook.edu for quick response to resolve any bug or feature update.**

### Installation
Please follow the following steps to install NanoBLASTer from source:
```
Clone NanoBLASTer source code: git clone https://github.com/ruhulsbu/NanoBLASTer.git
Go to the NanoBLASTer source directory: cd NanoBLASTer/nano_src 
Build the NanoBLASTer project: make
```  

### Input specifications
```
Use the following options to run NanoBLASTer:
-C: To specify one of the Parameters: -C10, -C25, or -C50
-r: To specify the name of Reference file (FASTA format)
-i: To specify the name of Reads file (FASTA format)
-o: To specify the prefix of Output file
-k: To specify the size of KMER
-a: To specify the size of ANCHOR
-l: To specify the min number of Clusters
-s: To run the program at higher sensitivity
-n: To specify the Number of reads to be aligned
-g: To specify the interval (or Gap) length between KMERs
-X: To configure NanoBLASTer for less memory using Single index
-h, or -?: To print this Help information.
```

### Usage examples
```
Run NanoBLASTer in ``fast'' mode (KMER=13, ANCHOR=45 and CLUSTERS=10):
$ ./nanoblaster -C10 -r path/to/reference.fa -i path/to/reads.fa -o output

Run NanoBLASTer in ``sensitive'' mode (KMER=11, ANCHOR=40 and CLUSTERS=25):
$ ./nanoblaster -C25 -r path/to/reference.fa -i path/to/reads.fa -o output

Run NanoBLASTer in ``highly sensitive'' mode (KMER=11, ANCHOR=40, CLUSTERS=50 and SENSITIVITY=TRUE):
$ ./nanoblaster -C50 -r path/to/reference.fa -i path/to/reads.fa -o output

Run NanoBLASTer with default KMER=11, ANCHOR=40 and CLUSTERS=10:
$ ./nanoblaster -r path/to/reference.fa -i path/to/reads.fa -o output

Run NanoBLASTer with KMER=13, ANCHOR=45 and CLUSTERS=25 using default parameters at higher sensitivity:
$ ./nanoblaster -r path/to/reference.fa -i path/to/reads.fa -o output -k13 -a45 -l25 -s

Run NanoBLASTer in ``fast'' mode using ``single index'':
$ ./nanoblaster -C10 -X -r path/to/reference.fa -i path/to/reads.fa -o output

Run NanoBLASTer in ``fast'' mode using ``interval=4'':
$ nanoblaster -C10 -g 4 -r path/to/reference.fa -i path/to/reads.fa -o output
```
* It is not recommended to use any additional parameters except input and output with C10, C25 or C50. But if we want to force additional parameters then we can add them after mentioning one of these sensitivity modes.

### Optimize configurations
Edit the configurations in constant.h file to optimize NanoBLASTer alignment manually. Editing the following constants will have an overall impact on the alignment quality of NanoBLASTer:
- To increase/decrease the gap penalty change the constant "GAP" (positive number, default 6)
- To increase/decrease the gap open penalty change the constant "GAP_OPEN" (additional penalty, default 2)
- To increase/decrease the base mismatch weight change the constant "MISMATCH" (negative number, default -4)
- To increase/decrease the base matching weight change the constant "WEIGHT" (positive number, default 5)
- To change the percentage of identity match for sequence alignment change the constant "KBAND_PERCENT_MATCH" (default 55)
- To change the size of read length considered for alignment change "MIN_READ_LEN" (default 100), or "MAX_READ_LEN" (default 200,000)
- If the genome is circular the constant "CIRCULAR" must be 0  

Editing the following constants will have an impact on the sensitivity and runtime of NanoBLASTer. These constants have been defined with the default values for NanoBLASTer. Some of these default values have been optimized for "fast", "sensitive" and "highly sensitive" mode. So these three configured modes will not be affected even if the following constants are edited. It is recommended to define -k, -a and -l when the following constants are changed manually:
- To increase/decrease the seed tuple distance in a cluster change the constant SEED_TUPLE_DIST
- To increase/decrease the primary anchor word size change the constant PRIM_ANCHOR_WORD
- To increase/decrease the secondary anchor word size change the constant SECND_ANCHOR_WORD
- To increase/decrease the primary anchor percent of identity change the constant PRIM_ANCHOR_PIDENT 
- To increase/decrease the secondary anchor percent of identity change the constant SECND_ANCHOR_PIDENT
- To increase/decrease the size of block for banded sequence alignment change the constant FRAGMENT_SIZE

Finally, define -k, -a and -l while running the ./nanoblaster executable with these changed constants to get the expected performance.

Currently, the sensitivity modes (C10, C25 and C25) are defined to align data set with very low accuracy (~65%). If the accuracy of the data set is ~85% then we need to modify the KMER and/or ANCHOR size when mentioning the sensitivity mode to get better runtime performance.

### Contact information
Please send your comments or bug reports to moamin@cs.stonybrook.edu 
