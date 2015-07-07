#include "library.h"

#ifndef STRUCTURE_H
#define STRUCTURE_H

typedef struct _node_ {
	int read_ind;
	int ref_ind;
	//int current_ind;
	//int weight;
	struct _node_ *next;
} node; 

typedef struct _reference_ {	//reference structure
	string ref;
	string rev;
	string name;
	//unordered_map<long, vector<int> > index;
	int *index;
	int *revind;
	vector<vector<int> > position;
} reference_index;

typedef struct _alignment_ {	//alignment structure
	vector<pair<char, char> > alignment;
	vector<pair<char, char> > end_to_end;
	vector<pair<int, int> > fragment_ind;
	int ref_start, ref_end;
	int read_start, read_end;
	int ref_ind, read_dir;
	int identity_match, total_len;
	int gaps, mismatches;
	int read_kmer, ref_kmer;
	int align_start, align_end;
} fragment_alignment;

typedef struct _cost_direction_ {
	int dir, cost;
	char ch1, ch2;
	int matrix_col;
	int str2_index;
	int match;
	int length;
} cell;


typedef vector<pair<char, char> > SeqVector;	//Sequence Vector
typedef pair<long long, SeqVector> SeqPair;	//Sequence Pair
typedef priority_queue<SeqPair> SeqPQ;		//Sequence Priority queue
/*
extern vector<reference_index> refindex;        //vector of reference
extern cell **matrix;// = NULL;                 //matrix for edit distance

extern int KMER;// = 11;                //Size of KMER. Choose all KMER from Chromosome and Target and hash them.
extern int SEED;// = 40;                //For the Set of KMER in due order find the SEED with 55% match
extern int MINREAD;// = 0;              //Start from Read No MINREAD
extern int MAXREAD;// = 100;            //End before MAXREAD is reached
extern int BOUNDARY;// = 10;            //Rank the alignment and take BOUNDARY option based on max SEED found
extern int MAX_MATCHED;// = 55;         //Maximum matching score
extern int MAX_SCORED;// = 0;

extern ofstream fp_error_free_seg;
extern ofstream fp_csv;
extern ofstream fp_blastn;
extern ostringstream logstr;

extern int t_lookup;// = 0;
extern int t_seed;// = 0;
extern int t_extend;// = 0;
extern int t_chain;// = 0;
extern int t_lis;// = 0;
extern int t_kmer_count;// = 0;
extern long base_power_kmer[12][4];
extern long error_free_seg[1000];
*/

#endif
