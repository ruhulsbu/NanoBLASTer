#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

void trace_path_cell_back(cell **path, int row, int column, string& str1, string& str2, vector<pair<char, char> >& alignment)
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

	//cout << "alignment validate/validation is done properly" << endl;
}

void print_alignment(vector<pair<char, char> >& alignment, string& ref, string& read, int ref_position, 
			int read_position, int step, fragment_alignment &fragment_alignment_info, bool print)
{
	int i, k, index = 0;
	int updated_index;
	int fragment_count = 0;
	int start = 0, end = -1;
	int  score = 0;
	bool flag_ref = false;
	bool flag_read = false;
	bool flag_break = false;
	int ref_start, ref_end;
	int read_start, read_end;

	ref_start = ref_position;//04-01-15
	read_start = read_position;

	while(index < alignment.size())
	{
		flag_break = false;
		if(print == true)
			cout << ref_position << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			/*
			if(alignment[i].second == '.')//03-26-15
			{
				updated_index = i;
				flag_break = true;
				break;
			}
			*/
			if(print == true)
				cout << alignment[i].first;
			
			if(alignment[i].first != '-')
			{
				/*
				if(flag_ref == false)
				{
					if(alignment[i].first != '-')//03-26-15
					//if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
					{
						start = i;
						flag_ref = true;
						ref_start = ref_position;
						ref_end = ref_position;
					}
				}
				else
				{
					if(alignment[i].first != '-')//03-26-15
					//if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
                                        {
						end = i;
                                                ref_end = ref_position;
                                        }

				}
				*/
				assert(alignment[i].first == ref.at(ref_position));//03-25-15
				ref_position += step;
			}
			
		}
		if(print == true)
		{
			//cout << "\t" << ref_position;
			printf("%12d", ref_position - step);
			cout << endl;

			cout << "\t    ";
		}
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			/*
			if(alignment[i].second == '.')
				break;
			*/
			if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
			{
				if(print == true)
					cout << "|";
				score += 1;
			}
			else
			{
				if(print == true)
					cout << " ";
			}
			
		}

		if(print == true)
		{
			cout << endl;
			cout << read_position << "\t    ";
		}
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			/*
			if(alignment[i].second == '.')//03-26-15
				break;
			if(alignment[i].second == 'o')
				cout << '-';
			else
				cout << alignment[i].second;
			*/
			if(print == true)
				cout << alignment[i].second;
			if(alignment[i].second != '-')
			{
				/*
				if(flag_read == false)
                                {
					if(alignment[i].first != '-')//03-26-15
                                        //if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
                                        {
                                                flag_read = true;
                                                read_start = read_position;
                                                read_end = read_position;
                                        }
                                }
                                else
                                {
					if(alignment[i].first != '-')//03-26-15
                                        //if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
                                        {
                                                read_end = read_position;
                                        }

                                }
				*/
				//cout << endl << alignment[i].second << " VS " << read.at(read_position) << endl;
				assert(alignment[i].second == read.at(read_position));
				read_position += step;
			}
			
		}

		if(print == true)
		{
			//cout << "\t" << read_position;
			
			printf("%12d", read_position - step);
			cout << endl;
			cout << endl;
		}

		if(flag_break == true)
		{
			index = updated_index;
			while(alignment[index].second == '.')
				index += 1;
			//while(ref_position < fragment_alignment_info.fragment_ind[fragment_count].first &&
			//	read_position < fragment_alignment_info.fragment_ind[fragment_count].second)
			fragment_count += 1;
			ref_position = fragment_alignment_info.fragment_ind[fragment_count].first;
			read_position = fragment_alignment_info.fragment_ind[fragment_count].second;
			cout << "Next ref_position = " << ref_position << ", And read_position = " << read_position << endl;
			cout << "BREAK##########################################################################DANCE" << endl << endl;
		}
		else
			index += BREAKAT;
	}

	ref_end = ref_position;//04-01-15
	read_end = read_position;
	start = 0;
	end = alignment.size() - 1;

	if(print == true)
	{
		cout << "##########################################################################" << endl;
		cout << "##########################################################################" << endl;
		cout << endl << endl;
		return;
	}

	fragment_alignment_info.ref_start = ref_start;
        fragment_alignment_info.ref_end = ref_end - 1;//<= end

	fragment_alignment_info.read_start = read_start;
	fragment_alignment_info.read_end = read_end - 1;//<= end

	for(i = start; i <= end; i++)
	{
		fragment_alignment_info.alignment.push_back(alignment[i]);
	}
	//fragment_alignment_info.end_to_end = alignment;
	fragment_alignment_info.identity_match = score;
	fragment_alignment_info.total_len = fragment_alignment_info.alignment.size();
}


