## NanoBLASTer: Basic Local Alignment and Search Tool for Oxford Nanopore Long Sequences  
**__Current Version: 0.16__**  
Release date: July 15, 2015  
Platform: Linux x64 system

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
-b: To specify the min number of Clusters (optional)
-n: To specify the Number of reads to be aligned (optional)
-h or -?: To print this Help information.
```

### Usage examples
```
# Align all reads from a given FASTA file to reference using given configuration:
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
- To increase/decrease the base mismatch weight change the constant "MISMATCH" (negative number)
- To increase/drease the base matching weight change the constant "WEIGHT"
- To increase/decrease the match open reward change the constant "MAT_OPEN"
- To change the size of band for block wise sequence alignment algorithm change the constant "FRAGMENT_SIZE"
- To change the percentage of identity match for sequence alignment change the constant "KBAND_PERCENT_MATCH"
- To change the size of read length considered for alignment change "MINREADLEN", or "MAXREADLEN"
- If the genome is circular the constant "BOUNDED" must be 0  

### Contact information
Please send your comments or bug reports to moamin@cs.stonybrook.edu 
