#include "library.h"
#include "constant.h"
#include "structure.h"

#ifndef GLOBAL_H
#define GLOBAL_H

extern int SINGLE;
extern int KMER;				//To specify the size of KMER 
extern int ANCHOR;				//To specify the ANCHOR size
extern int MINREAD;				//To specify the number of reads to be skipped
extern int MAXREAD;				//To specify maximum number of reads to be aligned
extern int CLUSTER;				//Maximum clusters to be used
extern int CLUSTERNUMBER;			//Clusters considered = CLUSTERNUMBER * CLUSTER;
extern int ALIGNMENT_CNT;			//Alignment Count for indels Analysis
extern bool HIGHSITIVE;				//To specify to run highly sensitive mode

//extern ofstream fp_error_free_seg;		//This file keeps the record of error free segment
extern ofstream fp_error_dist;			//This file keeps the record of error distribution
extern ofstream fp_csv;				//This csv file keep statistics for each read
extern ofstream fp_blastn;			//This files keeps the statistics for each alignment
extern ostringstream logstr;			//It stores the name of output

/*These variables are used keep the records of
time taken by each module of NanoBLASTer*/
extern int t_lookup;
extern int t_anchor;
extern int t_extend;
extern int t_chain;
extern int t_lis;
extern int t_kmer_count;

extern long base_power_kmer[21][4];		//Pre-caculated pow(BASE, i) * k for (i, k)
//extern long error_free_seg[1000];		//To calculate error free distribution
extern long error_dist[10];			//To calculate indel distribution

extern reference_index basic_index;
extern vector<reference_index> refindex;	//Vector of reference string
extern cell **matrix;// = NULL;			//Matrix for sequence similarity

extern int SEEDTUPLEDIST;
extern int MINREADLEN;
extern int MAXREADLEN;
extern int INTERVALX;

extern int PRIMANCHORWORD;
extern int SECNDANCHORWORD;

extern float PRIMANCHORPIDENT;
extern float SECNDANCHORPIDENT;

extern int FRAGMENTSIZE;
extern float ALIGNMENTCOVERAGE;
extern int ORIGINALPERCENTMATCH;
extern float EXPLICIT_SLOPE_90_CHEK;
extern float EXPLICIT_SLOPE_80_CHEK;
extern float EXPLICIT_SLOPE_70_CHEK;
extern float EXPLICIT_SLOPE_50_CHEK;
#endif
