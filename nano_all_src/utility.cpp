#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

int map_value(char ch)
{
	if(ch == 'A')
		return 0;
	if(ch == 'C')
		return 1;
	if(ch == 'G')
		return 2;
	if(ch == 'T')
		return 3;

	return -1;
}

void init_matrix()
{
        matrix = new cell *[2 * FRAGMENT_SIZE + 5];
        for(int i = 0; i < 2 * FRAGMENT_SIZE + 5; i++)
        {
                matrix[i] = new cell[FRAGMENT_SIZE + 5];
        }

}

void remove_matrix()
{
	for(int i = 0; i < 2 * FRAGMENT_SIZE + 5; i++)
        {
                delete [] matrix[i];
                matrix[i] = NULL;
        }

        delete [] matrix;

        matrix = NULL;
}

int similarity(char x, char y)
{
	if(x == y)
		return 1 * WEIGHT;
	else 
		return MISMATCH;
}

int min(int x, int y)
{
	if(x < y)
		return x;
	else 
		return y;
}

int max(int x, int y) 
{
	if(x < y)
		return y;
	else
		return x;
}

void upper_case(string& str)
{
	char *array = new char [str.length() + 1];
	strcpy(array, str.c_str());
	for(int i = 0; i < str.length(); i++)
		array[i] = toupper(array[i]);
	
	str.clear();
	str = string(array);
	delete [] array;
}

void reverse_str(string &str)
{
	int i, k;
	char ch;
	
	//cout << "Forward = " << str.substr(0, 80) << endl;
	for(i = 0, k = str.length() - 1; i < k; i++, k--)
	{
		ch = str[i];
		str[i] = str[k];
		str[k] = ch;
	}
	//cout << "Reverse = " << str.substr(str.length() - 80, 80) << endl;

}

string reverse_complement(string& str)
{
	int i, j;
	char x, y;
	char *input = new char [str.length() + 1];
	char *output = new char [str.length() + 1];
	strcpy(input, str.c_str());
	if(DEBUG == 1) cout << str << endl;	
	for(i = str.length() - 1, j = 0; i >= 0; i--, j++)
	{
		x = input[i];
		if(x == 'A')
			y = 'T';
		else if(x == 'T')
			y = 'A';
		else if(x == 'C')
			y = 'G';
		else if(x == 'G')
			y = 'C';
		else 
			y = 'N';//might also be = x here

		output[j] = y;
		if(DEBUG == 1) cout << x << " " << y << endl; 	
	}

	output[j] = '\0';
	if(DEBUG == 1) printf("%s\n", output);
	string rc(output);
	if(DEBUG == 1)
		cout << "Reverse(" << str << ") = " << rc << endl;
	
	delete [] input;
	delete [] output;

	return rc;
}

bool validate_alignment(vector<reference_index>& refindex, int ref_ind, vector<pair<char, char> >& alignment, 
			int ref_start, int read_start, string& read)
{
        for(int i = 0; i < alignment.size(); i++)
        {
                //cout << "Refs: " << ref_start << " = " << alignment[i].first << " VS " << refindex[ref_ind].ref.at(ref_start) << endl;
                //cout << "Read: " << read_start << " = " << alignment[i].second << " VS " << read.at(read_start) << endl;

                if(alignment[i].first != '-')
                {
			//cout << "Refs: " << ref_start << " = " << alignment[i].first << " VS " << refindex[ref_ind].ref.at(ref_start) << endl;
                        assert(alignment[i].first == refindex[ref_ind].ref.at(ref_start));
                        ref_start += 1;
                }
                if(alignment[i].second != '-')
                {
			//cout << "Read: " << read_start << " = " << alignment[i].second << " VS " << read.at(read_start) << endl;
                        assert(alignment[i].second == read.at(read_start));
                        read_start += 1;
                }
        }

	cout << "alignment validate/validation is done properly" << endl;
}

bool validate_chain(vector<reference_index>& refindex, int ref_ind, fragment_alignment& fragment_alignment_info, string& read)
{
	int ref_start = fragment_alignment_info.ref_start;
	int read_start = fragment_alignment_info.read_start;
	for(int i = 0; i < fragment_alignment_info.alignment.size(); i++)
	{
		//cout << ref_start << " = " << fragment_alignment_info.alignment[i].first << " VS " << refindex[ref_ind].ref.at(ref_start) << endl;
		//cout << read_start << " = " << fragment_alignment_info.alignment[i].second << " VS " << read.at(read_start) << endl;

		if(fragment_alignment_info.alignment[i].first != '-')
		{
			assert(fragment_alignment_info.alignment[i].first == refindex[ref_ind].ref.at(ref_start));
			ref_start += 1;
		}
		if(fragment_alignment_info.alignment[i].second != '-')
		{
			assert(fragment_alignment_info.alignment[i].second == read.at(read_start));
			read_start += 1;
		}
	}
	cout << "chain validate/validation is done properly" << endl;
}
