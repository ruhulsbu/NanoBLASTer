#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

void sam_format(fragment_alignment &final_alignment_info, vector<reference_index>& refindex, 
					string& read, string& read_name, vector<string>& output)
{
	//cout << "####################################################################################\n";
	//cout << "Final Result - " << ":\n";
	
	int read_dir = 0, ref_ind = 0, ref_position = 0, maximum = 0, editdist = 0;
	int match = 0, insert = 0, delet = 0, subst = 0, ignore = 0;
	int previousop = 0, currentop = 0;
	int total_len, total_score, ref_length;
	string cigar = "", cigarseq = "";

	vector<pair<char, char> >& alignment = final_alignment_info.alignment;

	if(!alignment.empty())
	{	
		ref_ind = final_alignment_info.ref_ind;
		read_dir = final_alignment_info.read_dir;
		total_score = final_alignment_info.identity_match;
		total_len = final_alignment_info.total_len;
		ref_position = final_alignment_info.ref_start; 
		ref_length = (final_alignment_info.ref_end - final_alignment_info.ref_start);

		//cout << "RefLen = " << refindex[ref_ind].ref.length() << ", ReadLen = " << read.length() << 
		//	" AlignmentLen = " << total_score << endl;

		ALIGNMENT_CNT += 1;
		memset(error_dist, 0, sizeof(long) * 10);

		error_dist[7] = total_len;
		error_dist[8] = read.length();
		error_dist[9] = ref_length;

		int total_match = 0;
		int total_insert = 0;
		int total_delete = 0;
		int total_subst = 0;
		int total_nchar = 0;
		int total_ignore = 0;
		int real_match = 0;
		
		for(int i = 0; i < alignment.size(); i++)
		{
			if(alignment[i].first == '-' && alignment[i].second == '.')
				continue;
			//if(alignment[i].first == '-' && alignment[i].second == 'N')
			//	continue;
			if(alignment[i].first == '-' && alignment[i].second == '-')
				continue;

			if(alignment[i].second == '.')
                        {
                                ignore += 1;
                                currentop = IGNORE;
                                total_ignore += 1;
                        }
			else if(alignment[i].first == alignment[i].second)
			{
				match += 1;
				real_match += 1;
				currentop = MATCH;
				total_match += 1;
				cigarseq += alignment[i].first;
			}
			else if(alignment[i].first != alignment[i].second && alignment[i].second == '-')
			{
				delet += 1;
				editdist += 1;
				currentop = DELETE;
				total_delete += 1;
			}
			else if(alignment[i]. first != alignment[i].second && alignment[i].first == '-')
			{
				insert += 1;
				editdist += 1;
				currentop = INSERT;
				total_insert += 1;
				cigarseq += alignment[i].second;
			}
			else if(alignment[i].first != alignment[i].second && alignment[i].second != '.')
			{	
				subst += 1;
				editdist += 1;
				currentop = SUBSTITUTE;
				if(alignment[i].second != 'N')
					total_subst += 1;
				else
					total_nchar += 1;
				cigarseq += alignment[i].second;
			}
			/*
			else if(alignment[i].first == '.' || alignment[i].second == '.')
			{
				ignore += 1;
				currentop = IGNORE;
				total_ignore += 1;
				cigarseq += alignment[i].second;
			}
			*/
			if(i == 0) previousop = currentop;

			if(currentop != previousop)
			{
				ostringstream numstr;
				if(previousop == MATCH)
				{
					numstr << match;
					cigar += numstr.str() + "M";
					//error_free_seg[real_match] += 1;
					error_dist[0] = max(error_dist[0], real_match);
					match = 0;
					real_match = 0;
				}
				else if(previousop == DELETE)
				{
					numstr << delet;
					cigar += numstr.str() + "D";
					delet = 0;
				}
				else if(previousop == INSERT)
				{
					numstr << insert;
					cigar += numstr.str() + "I";
					insert = 0;
				}
				else if(previousop == SUBSTITUTE)
				{
					numstr << subst;
					cigar += numstr.str() + "X";
					subst = 0; 
				}
				else if(previousop == IGNORE)
				{
					numstr << ignore;
					cigar += numstr.str() + "N";
					ignore = 0;
				}

				//cout << "Updating previos op = " << previousop << " with CIGAR = " << cigar << endl;
				previousop = currentop;

			}

			//cout << alignment[i].first << " and " << alignment[i].second << " where op = " << currentop << endl;
		}

		ostringstream numstr;
		if(previousop == MATCH)
		{
			numstr << match;
			cigar += numstr.str() + "M";
			//error_free_seg[real_match] += 1;
			error_dist[0] = max(error_dist[0], real_match);
			match = 0;
			real_match = 0;
		}
		else if(previousop == DELETE)
		{
			numstr << delet;
			cigar += numstr.str() + "D";
			delet = 0;
		}
		else if(previousop == INSERT)
		{
			numstr << insert;
			cigar += numstr.str() + "I";
			insert = 0;
		}
		else if(previousop == SUBSTITUTE)
		{
			numstr << subst;
			cigar += numstr.str() + "X";
			subst = 0;
		}
		else if(previousop == IGNORE)
		{
			numstr << ignore;
			cigar += numstr.str() + "N";
			ignore = 0;
		}


		error_dist[1] = total_match;
		error_dist[2] = total_insert;
		error_dist[3] = total_delete;
		error_dist[4] = total_subst;
		error_dist[5] = total_nchar;
		error_dist[6] = total_ignore;

		if(DEBUG == 99)
		{
			fp_error_dist << ALIGNMENT_CNT << "," << error_dist[0] << "," << error_dist[1] <<
					"," << error_dist[2] << "," << error_dist[3] <<
					"," << error_dist[4] << "," << error_dist[5] <<
					"," << error_dist[6] << "," << error_dist[7] <<
					"," << error_dist[8] << "," << error_dist[9] << endl;
					//max_error_seg, total_match, total_insert, total_delete
					//total_substitute, total_nchar, total_ignore, alignment_length,
					//read_length, reference_length
		}
		/*
		cout << "total_match = " << total_match << endl;
		cout << "total_insert = " << total_insert << endl;
		cout << "total_delete = " << total_delete << endl;
		cout << "total_subst = " << total_subst << endl;
		cout << "total_nchar = " << total_nchar << endl;
		cout << "total_ignore = " << total_ignore << endl;
		cout << "total reference = " << ref_length << endl;
		*/
		//assert(cigarseq.length() == (alignment.size() - total_delete));
		assert(total_match + total_insert + total_delete + total_subst + total_nchar + total_ignore == 
					total_insert + ref_length);

		//cout << "Chain length = " << alignment.size()  << ", and Edit Distance: " << editdist <<  
		//", In Direction: " << read_dir  << " at ref_positin = " << ref_position << ", while CIGAR = " << cigar.length() << endl;
		if(DEBUG == 99)
		{
			fp_csv << read_dir << ", " << refindex[ref_ind].name << ", " << refindex[ref_ind].ref.length() << ", ";
			fp_csv << ref_position << ", " << total_score << ", " << total_len << ", ";
			fp_csv << ((100.00 * total_score) / (1.00 * total_len)) << ", " << alignment.size() << ", ";
			fp_csv << (1.0 * total_len / read.length()) << ", ";	
		}
	}
	else
	{
		total_score = alignment.size();	
		fp_csv << "0, 0, 0" << ", " << ref_position << ", " << total_score << ", 0, 0, 0, 0, ";
		//cout << "####################################################################################\n";
		return;	
	}

	output.push_back(read_name);			//1. read_name
	if(total_score == 0)				//2. direction
		output.push_back("4");
	else if(read_dir == 1)
		output.push_back("0");
	else
		output.push_back("16");
	output.push_back(refindex[ref_ind].name);	//3. reference name
	ostringstream numstr;				//4. index
	numstr << (ref_position + 1);		
	output.push_back(numstr.str());
	output.push_back("255");			//5. default mapq previously *
	if(cigar.length() == 0)				//6. cigar
		output.push_back("*");	
	else
		output.push_back(cigar);
	output.push_back("*");				//7. rnext
	output.push_back("0");				//8. pnext
	output.push_back("0");				//9. tlen 
	if(cigarseq.length() != 0)			//10. seq as cigar
		output.push_back(cigarseq);
	else 
		output.push_back("*");
	output.push_back("*");				//11. default quality

	//cout << "####################################################################################\n";
	//cout << endl << endl;

	return;
}

