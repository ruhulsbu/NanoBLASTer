#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

int create_gap_alignment(vector<pair<char, char> >& alignment)
{
	int mismatch = 0;
	int score = 0;
	int SEEDS = 100;
	int i;
	/*
	for(int i = 0; i < alignment.size() - SEEDS + 1; i++)
	{
		mismatch = 0;
		for(int k = 0; k < SEEDS; k++)
		{
			if(alignment[i + k].first != alignment[k + i].second)
				mismatch += 1;
		}

		if(1.0 * mismatch / (1.0 * SEEDS) < 0.55)
		{
			for(int k = 0; k < SEEDS; k++)
	                {
        	        	alignment[i + k] = make_pair(alignment[k + i].first, '-');

                	}
		}
		else
		{
			for(int k = 0; k < SEEDS; k++)
		        {
                	        if(alignment[i + k].first == alignment[k + i].second)
                        	        score += WEIGHT;
                	}
		}

		i = i + SEEDS - 1;
	}
	*/
	/*
	int i, k, count;
	for(i = alignment.size() - 1; i > alignment.size() / 2 && i > WINDOW; i--)
	{
		count = 0;
		for(k = 0; k < WINDOW; k++)
		{
			if(alignment[i - k].first != alignment[i - k].second)
				break;
			count += 1;
		}
		if(count == WINDOW)
			break;
	}

	if(i + 1 > 0 && i + 1 < alignment.size())
                alignment.erase(alignment.begin() + i + 1, alignment.end());
	*/

	/*
	//useful
	int i, k, count;
	cout << "size of KMER = " << KMER << endl;
        for(i = 0; i < alignment.size() / 1.5 - WINDOW; i++)
        {
                count = 0;
                for(k = 0; k < KMER; k++)
                {
                        if(alignment[i + k].first != alignment[i + k].second)
                                break;
                        count += 1;
                }
                if(count == KMER)
                        break;
        }

	if(i > 0 && i < alignment.size())
                alignment.erase(alignment.begin(), alignment.begin() + i);
	*/

	for(i = 0; i < alignment.size(); i++)
        {
                //cout << alignment[i].first << ", " << alignment[i].second << endl;

                if(alignment[i].first == alignment[i].second)
                {
                        if(alignment[i].first != '-')
                                score += 1;
                }
	}

	cout << "unoptimized length = " << alignment.size() << ", and score = " << score << endl;

	if(100.0 * score / alignment.size() > 1.0 * GAP_PERCENT_MATCH)
		return alignment.size();

	for(i = 0; i < alignment.size(); i++)
	{
		if(alignment[i].first == alignment[i].second)
		{
			if(alignment[i].first != '-')
			{
				if(100.0 * score / (alignment.size() - i - 1) > 1.0 * GAP_PERCENT_MATCH)
					break;

				score = score - 1;
			}

		}
	}	

	if(i > 0 && i <= alignment.size())
	{
                //alignment.erase(alignment.begin(), alignment.begin() + i);
		for(int k = 0; k < i; k++)
		{
			if(alignment[k].second != '-')
				alignment[k] = make_pair(alignment[k].first, 'o');
			else
				alignment[k] = make_pair(alignment[k].first, '-');
		}
	}
	
	cout << "alignment length = " << alignment.size() << ", optimized length = " << (alignment.size() - i) << 
			", and score = " << score << endl;

	return (alignment.size() - i);
}

int optimize_path(vector<pair<char, char> >& alignment)
{
	int match = 0, i;
	for(i = 0; i < alignment.size(); i++)
	{
		//cout << alignment[i].first << ", " << alignment[i].second << endl;

		if(alignment[i].first == alignment[i].second) 
		{
			if(alignment[i].first != '-')
				break;
			else
				match = match + MISMATCH;
		}
		else
			match = match - GAP;

	}

	if(i > 0 && i < alignment.size())
		alignment.erase(alignment.begin(), alignment.begin() + i);
	
	for(i = alignment.size() - 1; i >= 0; i--)
	{
		//cout << alignment[i].first << ", " << alignment[i].second << endl;

		if(alignment[i].first == alignment[i].second) 
		{
			if(alignment[i].first != '-')
				break;
			else
				match = match + MISMATCH;
		}
		else
			match = match - GAP;

	}

	if(i + 1 > 0 && i + 1 < alignment.size())
		alignment.erase(alignment.begin() + i + 1, alignment.end());

	for(i = 0; i < alignment.size(); i++)
	{
		//cout << alignment[i].first << ", " << alignment[i].second << endl;

		if(alignment[i].first == alignment[i].second) 
		{
			if(alignment[i].first != '-')
				match = match + WEIGHT;
			else
				match = match + MISMATCH;
		}
		else
			match = match - GAP;

	}


	return match;
}

void print_path_matrix(long long **path, int row, int column)
{
	int i, j;
	cout << "ROW = " << row << " and COL = " << column << endl;
	if(path == NULL)
	{
		cout << "NULL Pointer Detected" << endl;
		return;
	}
	for(i = 0; i < row; i++)
	{
		for(j = 0; j < column; j++)
		{
			cout << "(" << i << "," << j << "): " << path[i][j] << "\t";
		}
		cout << endl;
	}

}

void print_path_cell_back(cell **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment)
{
	int dir, k;
	char x, y;

	while(row != 0 || column != 0)
	{	
		dir = (int) path[row][column].dir;
		k = path[row][column].str2_index;

		if(dir == DIAG)
		{
			x = str1.at(row - 1);
			y = str2.at(k - 1);
		}
		else if(dir == UP)
		{
			x = str1.at(row - 1);
			y = '-';
		}
		else
		{
			x = '-';
			y = str2.at(k - 1);
		}

		//if(!(x == '-' && y == 'N'))
		alignment.push_back(make_pair(x, y));


		column = path[row][column].matrix_col;

		if(dir == DIAG || dir == UP) 
			row = row - 1;	
	}

}



void print_path_back(long long **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment)
{
	int dir, k;
	char x, y;

	while(row != 0 || column != 0)
	{	
		dir = (int) ((path[row][column] / MAXLEN) / MAXLEN);
		k = path[row][column] % MAXLEN;

		if(dir == DIAG)
		{
			x = str1.at(row - 1);
			y = str2.at(k - 1);
		}
		else if(dir == UP)
		{
			x = str1.at(row - 1);
			y = '-';
		}
		else
		{
			x = '-';
			y = str2.at(k - 1);
		}

		//if(!(x == '-' && y == 'N'))
		alignment.push_back(make_pair(x, y));


		column = (path[row][column] / MAXLEN) % MAXLEN;

		if(dir == DIAG || dir == UP) 
			row = row - 1;	
	}

}

void print_path(long long **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment)
{
	int i, j, k, dir;
	char x, y;
	
	i = row;
	k = path[row][column] % MAXLEN;
	j = (path[row][column] / MAXLEN) % MAXLEN;
	dir = (int) ((path[row][column] / MAXLEN) / MAXLEN);


	assert(i >= 0);
	assert(k >= 0);
	assert(j >= 0);
	assert(dir >= 0);
	assert(i <= str1.length());
	assert(k <= str2.length());
	assert(j <= 2 * KBAND + 2);
	
	if(DEBUG == 1) cout << row << "," << column << ": " << path[row][column] << endl;
	if(row == 0 && column == 0)
		return;
	
	if(DEBUG == 1)
		cout << "i = " << i << ", k = " << k << ", j = " << j << endl;

	//if(i == 0 || k == 0)
	//	return;

	if(dir == DIAG)
		print_path(path, i - 1, j, str1, str2, alignment);
	else if(dir == UP)
		print_path(path, i - 1, j, str1, str2, alignment);
	else
		print_path(path, i, j, str1, str2, alignment);
	
	if(dir == DIAG)
	{
		//if(str1.at(i - 1) == str2.at(k - 1))
		//	x = y = str1.at(i - 1);

		x = str1.at(i - 1);
		y = str2.at(k - 1);
	}
	else if(dir == UP)
	{
		x = str1.at(i - 1);
		y = '-';
	}
	else
	{
		x = '-';
		y = str2.at(k - 1);
	}

	alignment.push_back(make_pair(x, y));
}


