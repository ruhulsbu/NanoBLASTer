#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

void print_alignments(vector<pair<char, char> >& alignment)
{
	int i, k, index = 0;

	while(index < alignment.size())
	{
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			cout << alignment[i].first;
		}

		cout << endl;

		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
				cout << "|";
			else
				cout << " ";
		}

		cout << endl;

		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			cout << alignment[i].second;
		}

		cout << endl << endl;
		index += BREAKAT;
	}
}

void print_alignment_back(vector<pair<char, char> >& alignment_given, int ref_position, int read_position, int step)
{
	int i, k, index = 0;
	pair<char, char> swap;
	vector<pair<char, char> > alignment;	

	for(k = alignment_given.size() - 1; k >= 0; k--)
	{
		alignment.push_back(alignment_given[k]);
	}
	
	while(index < alignment.size())
	{
		cout << ref_position << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			cout << alignment[i].first;
			if(alignment[i].first != '-')
				ref_position += step;
		}
		//cout << "\t" << ref_position;
		printf("%12d", ref_position - step);
		cout << endl;

		cout << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
				cout << "|";
			else
				cout << " ";
		}

		cout << endl;

		cout << read_position << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].second == 'o')
				cout << '-';
			else
				cout << alignment[i].second;
			if(alignment[i].second != '-')
				read_position += step;
		}
		//cout << "\t" << read_position;
		printf("%12d", read_position - step);
		cout << endl;
		cout << endl;
		index += BREAKAT;
	}
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
		cout << ref_position << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].second == '.')//03-26-15
			{
				updated_index = i;
				flag_break = true;
				break;
			}
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
		//cout << "\t" << ref_position;
		printf("%12d", ref_position - step);
		cout << endl;

		cout << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].second == '.')
				break;
			if(alignment[i].first == alignment[i].second && alignment[i].first != '-')
			{
				cout << "|";
				score += 1;
			}
			else
				cout << " ";
		}

		cout << endl;

		cout << read_position << "\t    ";
		for(i = index, k = 0; i < alignment.size() && k < BREAKAT; i++, k++)
		{
			if(alignment[i].second == '.')//03-26-15
				break;
			if(alignment[i].second == 'o')
				cout << '-';
			else
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
		//cout << "\t" << read_position;
		printf("%12d", read_position - step);
		cout << endl;
		cout << endl;
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
		return;

	fragment_alignment_info.ref_start = ref_start;
        fragment_alignment_info.ref_end = ref_end - 1;//<= end

	fragment_alignment_info.read_start = read_start;
	fragment_alignment_info.read_end = read_end - 1;//<= end

	for(i = start; i <= end; i++)
	{
		fragment_alignment_info.alignment.push_back(alignment[i]);
	}
	fragment_alignment_info.end_to_end = alignment;
	fragment_alignment_info.identity_match = score;
	fragment_alignment_info.total_len = fragment_alignment_info.alignment.size();
}


