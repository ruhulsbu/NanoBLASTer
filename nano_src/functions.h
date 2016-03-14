#include "library.h"
#include "constant.h"
#include "structure.h"

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

//main.cpp
void input_reference(vector<pair<string, string> >& reference, string& ref_file, string& sam_file);
//end

//align_reads.cpp
extern void align_reads(vector<pair<string, string> >& reference, string& read_file, string& sam_file, vector<reference_index>& refindex);
extern void read_vs_reference(string& read, string& read_name, int dir, vector<reference_index>& refindex, 
			vector<pair<int, vector<pair<int, int> > > >& primary_chain);
extern void refine_kmer_index(vector<pair<int, int> >& kmer_index, vector<pair<int, vector<pair<int, int> > > >& primary_chain,
                                string read, int dir, vector<reference_index>& refindex, int index);
extern void create_primary_chain_from_list(node *current_node, vector<pair<int, int> >& chain, node *next_node, node *previous_node,  int readlen);
//end

//alignment.cpp
extern void align(string& read, string& read_name, int direction, vector<reference_index>& refindex, 
			vector<pair<int, vector<pair<int, int> > > >& primary_chain, vector<vector<string> >& output);
extern void create_alignment(string& read, string& reverse, string& read_name, vector<reference_index>& refindex, 
	vector<pair<int, pair<int, int> > >& kmer_ref, vector<fragment_alignment>& list_alignment_info);
extern void create_chain_alignment(string& read, int read_dir, vector<reference_index>& refindex, int ref_ind, vector<pair<int, int> >& chain,
                vector<fragment_alignment>& fragment_chain, int &total_len, int &total_score);
extern void create_final_alignment(string& read, int read_dir, vector<reference_index>& refindex, int ref_ind, vector<pair<int, int> >& chain,
                vector<fragment_alignment>& fragment_chain, int &total_len, int &total_score);
//end

//seed_selection.cpp
extern void find_lis_vector(vector<int>& lis_vector, vector<int>& indexing);
extern int count_kmer_match(string& str1, string& str2, int start, int end);
extern int find_hband_similarity(string& str1, string& str2, priority_queue<long>& pqueue, int start, int end);
//end

//edit_distance
extern int find_kband_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, int &str1_end, int &str2_end);
extern int find_banded_similarity(string& str1, string& str2, vector<pair<char, char> >& alignment, 
					int &str1_end, int &str2_end, bool direction);
//end

//optimize_kband.cpp
extern void trace_path_cell_back(cell **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment);
extern bool validate_alignment(vector<reference_index>& refindex, int ref_ind, vector<pair<char, char> >& alignment, 
				int ref_start, int read_start, string& read);
extern void print_alignment(vector<pair<char, char> >& alignment, string& ref, string& read, int ref_position, 
		int read_position, int step, fragment_alignment &fragment_alignment_info, bool flag);
//end

//sam_format.cpp
extern void sam_format(fragment_alignment &final_alignment_info, vector<reference_index> &refindex, string& read, 
		string& read_name, vector<string>& output);
//end

//utility.cpp
extern int map_value(char ch);
extern void init_matrix();
extern void remove_matrix();

extern int similarity(char x, char y);
extern int min(int x, int y);
extern int max(int x, int y);

extern void upper_case(string& str);
extern void reverse_str(string &str);
extern string reverse_complement(string& str);
extern int getValue();
//end

#endif
