#include "library.h"
#include "constant.h"
#include "structure.h"

#ifndef GLOBAL_H
#define GLOBAL_H

extern int KMER;// = 11;		//Size of KMER. Choose all KMER from Chromosome and Target and hash them. 
extern int SEED;// = 40;		//For the Set of KMER in due order find the SEED with 55% match
extern int MINREAD;// = 0;		//Start from Read No MINREAD
extern int MAXREAD;// = 100;		//End before MAXREAD is reached
extern int BOUNDARY;// = 10;		//Rank the alignment and take BOUNDARY option based on max SEED found
extern int MAX_MATCHED;// = 55;		//Maximum matching score
extern int MAX_SCORED;// = 0;
extern int ALIGNMENT_CNT;// = 0;	//Alignment Count for Insert, Delete, Substitute Analysis
extern bool HIGHSITIVE;

extern ofstream fp_error_free_seg;
extern ofstream fp_error_dist;
extern ofstream fp_csv;
extern ofstream fp_blastn;
extern ostringstream logstr;

extern int t_lookup;// = 0;
extern int t_seed;// = 0;
extern int t_extend;// = 0;
extern int t_chain;// = 0;
extern int t_lis;// = 0;
extern int t_kmer_count;// = 0;
extern long base_power_kmer[21][4];
extern long error_free_seg[1000];
extern long error_dist[10];

extern vector<reference_index> refindex;	//vector of reference
extern cell **matrix;// = NULL;			//matrix for edit distance

/*
extern int map_value(char ch);
extern void init_matrix();
extern void remove_matrix();
extern int similarity(char x, char y);
extern int min(int x, int y);
extern int max(int x, int y);

extern int is_overlap(string& str1, string& str2, SeqPQ& sequences, int dir);
extern int create_gap_alignment(vector<pair<char, char> >& alignment);

extern int optimize_path(vector<pair<char, char> >& alignment);
extern void print_path_matrix(long long **path, int row, int column);
extern void print_path_back(long long **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment);
extern void print_path(long long **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment);

extern int find_kband_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int &str1_end, int &str2_end);
extern int find_banded_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int &str1_end, int &str2_end);
extern int find_hband_similarity(string& str1, string& str2, priority_queue<long>& pqueue, int start, int end);

extern void print_alignment_back(vector<pair<char, char> >& alignment, int ref_position, int read_position, int step);
extern void print_alignment(vector<pair<char, char> >& alignment, int ref_position, int read_position,
                int step, fragment_alignment &fragment_alignment_info, bool flag);
*/
#endif
