## NanoBLASTer: Basic Local Alignment and Search Tool for Oxford Nanopore Long Sequences  
**__Current Version: 0.16__**  
Release date: July 15, 2015  
Platform: Linux x64 system

### Introduction
The quality of the Oxford Nanopore MinION sequenced DNA reads has been considered as the limiting factor for the success of nanopore sequencing technology. Thus, it has been of great interest to the researchers to build software tools that can analyze the long noisy reads at higher sensitivity. So far, the existing alignment methods like LAST and BLAST have been used to analyze these sequences, but these tools are unable to produce long alignments for noisy reads, and take too much processing time at higher sensitivity. To address these issues, we introduce a novel alignment tool, NanoBLASTer for the long noisy nanopore reads. This tool produces long alignments of the noisy nanopore reads at higher sensitivity. The run time of NanoBLASTer remains feasible even as we increase the sensitivity of the tool. 

### Methods
NanoBLASTer uses fixed size exact matching seeds followed by DP-based extension. However, because of the high error rate of the nanopore sequencing instruments (approximately 10% to 50% base error rate), the seeds that must be used are extremely short and provide relatively little specificity. NanoBLASTer overcomes this challenge and improves the specificity of short seeds by clustering neighboring seeds into mapping regions, and then identifying highly similar segment that we call ANCHORs from the clustered seeds. Extending the top scoring candidate ANCHORs with a block-wise banded sequence alignment algorithm generates the alignments. NanoBLASTer aligns long noisy reads using the following steps:
- Create an inverted K-mer index of the reference genome
- Identify and cluster neighboring K-mers seeds
- Identify and score the candidate ANCHOR
- Extend the ANCHOR into a complete alignment
- Report the alignments of sufficient quality in SAM format

### Results
- Aligning nanopore reads accurately has proven difficult because of the high base error rates. To mitigate the problem,
we present a novel long noisy read aligner, NanoBLASTer, custom designed for MinION sequenced reads. This tool can
align noisy reads at higher sensitivity. In our experiments of aligning noisy yeast reads sequenced with MinION, we found
that NanoBLASTer can generate alignments that are more than one kilobase longer on average than LAST.
- The run-time performance of NanoBLASTer does not degrade as we increase the sensitivity of the tool. This is not true in
case of LAST or BLAST. In our experiment, we run LAST and NanoBLASTer in their most sensitive configuration to align
the reads of yeast sequenced using MinION. NanoBLASTer runs at least 13 times faster than LAST (in most sensitive
configuration).
- Concerns have been raised about the base-level accuracy of nanopore sequences, but quantifying its performance is difficult
in the absence of a sensitive alignment program. To shed insights into the accuracy of the sequencer, we analyze several
characteristics of noisy sequences. For a particular yeast data set, we show that approximately 25% of the sequences align
with an average alignment length of 6:45kbp at 62:3% accuracy providing an average coverage depth of 11.72. We also
show that, for the long reads (> 40kb), on average half of their length can be aligned. Thus the long reads of MinION can be
used to span repetitive elements of sequences to provide highly contiguous assemblies.

**Preprint**: [NanoBLASTer: Fast Alignment and Characterization of Oxford Nanopore Single Molecule Sequence Reads]()  

**Nanopore sequencing data** of Yeast: [Yet to be published]()  
  
### Installation
Please follow the following steps to install NanoBLASTer from source:
```
Download NanoBLASTer source code: wget https://github.com/ruhulsbu/LongReadAlignment/archive/master.zip
Unzip NanoBLASTer project: unzip master.zip
Go to the NanoBLASTer directory: cd Project_directory 
Build the NanoBLASTer project: make
```  

### Input specifications
```
Use the following options to run NanoBLASTer:
-p: To specify one of the Configurations: C10, C25, or C50
-r: To specify the name of Reference file in FASTA format
-i: To specify the name of Reads file in FASTA format
-o: To specify the name of Output file
-k: To specify the size of KMER (optional)
-s: To specify the size of SEED (optional)
-b: To specify the min number of Trials (optional)
-n: To specify the Number of reads to be aligned (optional)
-h or -?: To print this Help information.
```

### Usage examples
```
# Align all reads from a given FASTA file to reference using given configurations:
time ./nano -r ../ref.fa -i ../reads.fa -p C10 -o output
- Gives all the aligned reads in sam format in output.sam file

# Align all reads from a given FASTA file to reference by user specified configurations:  
time ./nano -r ../ref.fa -i ../reads.fa -k 11 -s 40 -b 25 -o output
- Aligning using KMER=11, SEED=40 and CLUSTERS=25
```

### Optimize configurations
Edit the configurations in constant.h file to optimize NanoBLASTer alignment manually:
- To increase/decrease the gap penalty change the constant "GAP" (positive number)
- To increase/decrease the gap open penalty change the constant "GAP_OPEN"
- To increase/decrease the base mismatch weight change the constant "MISMATCH"
- To increase/drease the base matching weight change the constant "WEIGHT"
- To increase/decrease the match open reward change the constant "MAT_OPEN"
- To change the size of band for block wise sequence alignment algorithm change the constant "FRAGMENT_SIZE"
- To change the percentage of identity match for sequence alignment change the constant "KBAND_PERCENT_MATCH"
- To change the size of read length considered for alignment change "MINREADLEN", or "MAXREADLEN"
- If the genome is circular the constant "BOUNDED" must be 0  

### Contact information
Please send your comments or bug reports to moamin@cs.stonybrook.edu 
