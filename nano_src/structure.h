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

typedef struct _reference_ {			//reference index
	string ref;
	string rev;
	string name;
	//unordered_map<long, vector<int> > index;
	int *index;
	int *revind;
	vector<vector<int> > position;
} reference_index;

typedef struct _alignment_ {			//alignment information
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

typedef struct _cost_direction_ {		//matrix cell information
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

#endif
