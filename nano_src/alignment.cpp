#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

bool compare_fragments(fragment_alignment first, fragment_alignment second)
{
	if(first.ref_start != second.ref_start)
        	return (first.ref_start < second.ref_start);
	/*	
	if(first.read_start != second.read_start)
		return (first.read_start < second.read_start);
	*/
	if(first.ref_end != second.ref_end)
		return (first.ref_end > second.ref_end);
	/*
	if(first.read_end != second.read_end)
		return (first.read_end < second.read_end);
	*/
	return (first.alignment.size() > second.alignment.size());
}

bool compare_chain(pair<int, vector<pair<int, int> > > first_pair, pair<int, vector<pair<int, int> > > second_pair)
{
        return (first_pair.second.size() > second_pair.second.size() );
}

void align(string& read, string& read_name, int direction, vector<reference_index>& refindex, 
	vector<pair<int, vector<pair<int, int> > > >& primary_chain, vector<vector<string> >& final_result)	
{
	int ref_alignment_count[refindex.size()];
	memset(ref_alignment_count, 0, sizeof(int) * refindex.size());

	string reverse = reverse_complement(read);

	sort(primary_chain.begin(), primary_chain.end(), compare_chain);
	/*
	for(int i = 0; i < primary_chain.size(); i++)
        {
                cout << "\t" << i << ") " << primary_chain[i].first << " == " << primary_chain[i].second.size() << endl;
                for(int k = 0; k < primary_chain[i].second.size(); k++)
                        cout << "\t\tref = " << primary_chain[i].second[k].first << " : read = " <<
                                primary_chain[i].second[k].second << " : diff = "  <<
                                (primary_chain[i].second[k].first - primary_chain[i].second[k].second) << endl;
        }
	*/

	int loop_limit = CLUSTER;
	//int loop_limit = refindex.size() * CLUSTER;	
	int first_index, last_index, index_dir, index;
	float matching_percentage, coverage_percentage, final_coverage = 0.0;
	vector<fragment_alignment> list_alignment_info;
	fragment_alignment final_alignment_info;
	
	for(int i = 0; i < primary_chain.size() && i < loop_limit; i++)
	{
		//if(i > 0 && primary_chain[i].second.size() * 2 < primary_chain[0].second.size())
		//	break;

		vector<pair<int, pair<int, int> > > kmer_ref;
		index_dir = primary_chain[i].first;
		index = index_dir / MULTIPLIER;
		//cout << "Reference Index = " << index << endl;
		
		vector<pair<int, int> > kmer_index = primary_chain[i].second;
		//cout << "Size of Kmer_Index = " << kmer_index.size() << endl;	

		for(int k = 0; k < kmer_index.size(); k++)
                {
                        if(kmer_index[k].first == -1)
                                continue;
			
			if(k > 0 && kmer_index[k].first - kmer_index[k - 1].first == INTERVALX &&
                                kmer_index[k].second - kmer_index[k - 1].second == INTERVALX)
                                continue;
			
                        first_index = kmer_index[k].first + 1 - kmer_index[k].second;///1.3
                        last_index = kmer_index[k].first + 1 - kmer_index[k].second + read.length() / 1.35;

                        //if((first_index >= 1) && (last_index <= refindex[index].ref.length()))
                        {
                                kmer_ref.push_back(make_pair(index_dir, kmer_index[k]));
                                //cout << "Pushing: first = " << kmer_index[k].first << " And second = " << kmer_index[k].second << endl;
                        }

		}
		
		if(kmer_ref.size() < 1)
		{
			kmer_index.clear();
			continue;
		}
		
		create_alignment(read, reverse, read_name, refindex, kmer_ref, list_alignment_info);

		if(list_alignment_info.size() > 0)
			loop_limit = CLUSTERNUMBER * CLUSTER - 2;

		kmer_ref.clear();
		kmer_index.clear();

	}

	for(int i = 0; i < list_alignment_info.size(); i++)
	{
		if(DEBUG == 99)
		{	
			final_alignment_info = list_alignment_info[i];
			fp_blastn << read_name << "," << refindex[final_alignment_info.ref_ind].name << ",";
			fp_blastn << KMER << "," << final_alignment_info.read_start << "," << final_alignment_info.read_end << ",";
			fp_blastn << final_alignment_info.ref_start << "," << final_alignment_info.ref_end << ",";
			fp_blastn << final_alignment_info.gaps << "," << final_alignment_info.mismatches << ",";
			fp_blastn << final_alignment_info.identity_match << "," << final_alignment_info.total_len << ",";
			fp_blastn << (100.00 * final_alignment_info.identity_match / final_alignment_info.total_len) << ",";
			fp_blastn << final_alignment_info.read_dir << "," << read.length() << ",";
			fp_blastn << (100.00 * final_alignment_info.total_len / read.length()) << endl;
		}
		vector<string> sam_output;
		//sam_format(final_alignment_info, refindex, read, read_name, sam_output);
		sam_format(list_alignment_info[i], refindex, read, read_name, sam_output);
		final_result.push_back(sam_output);
		//break;
	}
}

void create_alignment(string& read, string& reverse, string& read_name, vector<reference_index>& refindex, 
	vector<pair<int, pair<int, int> > >& kmer_ref, vector<fragment_alignment>& list_alignment_info)
{	
	//cout << "\tIn function align " << read_name << ": has length = " <<  read.length() << endl;

	long newentry;
	int ref_ind, read_dir;
	int index, direction, start, end = 0;
	int ret_val, pre_str_pos, pre_ref_pos;

	deque<pair<long, long> > match_queue;
	vector<pair<int, int> > chain;

	//clock_t t_start, t_end; 
	//time_t tstart, tend;

	//time(&tstart);
	//t_start = clock();
	string xstring, ystring;
	string strfrd, strrev;

	index = kmer_ref[0].first / MULTIPLIER;
	direction = kmer_ref[0].first % MULTIPLIER;

	for(int k = 0; k < kmer_ref.size(); k++)
	{
			
		if(k > 0)//03-26-15
		{
			if(kmer_ref[k].second.first - end < FRAGMENTSIZE)
				continue;
		}
		
		start = kmer_ref[k].second.first + 1;
		end = kmer_ref[k].second.first + 1;//second
		priority_queue<long> pqueue;
		/*		
		cout << "\t" << k << ")comparing with ref_ind = " << index << ", in the direction " << direction
			<< ", at starting positin = " << start << ", and end = " << (end + ANCHOR) << endl;

		cout << "reference length = " << refindex[index].ref.length() << ", and read length = " << read.length() 
			<< " from index = " << kmer_ref[k].second.second << endl;
		*/
		ret_val = -1;

		if(direction == FF && start + ANCHOR <= refindex[index].ref.length())
		{
			strfrd = read.substr(kmer_ref[k].second.second, min(ANCHOR, read.length() - kmer_ref[k].second.second));
			ret_val = find_hband_similarity(refindex[index].ref, strfrd, pqueue, start,
					min(refindex[index].ref.length(), end + ANCHOR));//read, max(1, start - ANCHOR)
		
			if(ret_val == -1)
			{
				int start1 = kmer_ref[k].second.first;
				int start2 = kmer_ref[k].second.second;
				int kmer_match = 0;
				while(refindex[index].ref.at(start1) == read.at(start2))
				{
					kmer_match += 1;
					start1 += 1;
					start2 += 1;
					if(start1 >= refindex[index].ref.length() || start2 >= read.length())
						break;
				}
				pre_str_pos = kmer_ref[k].second.second + kmer_match - ANCHOR;
				pre_ref_pos = start + kmer_match - ANCHOR;
				if(pre_str_pos >= 0 && pre_ref_pos >= 1)
				{		
					//cout << "max total number of character could be matched = " << kmer_match << endl;
					strfrd = read.substr(pre_str_pos, min(ANCHOR, read.length() - pre_str_pos));
					ret_val = find_hband_similarity(refindex[index].ref, strfrd, pqueue, pre_ref_pos,
                                		min(refindex[index].ref.length(), pre_ref_pos + ANCHOR));
					//cout << "backward = " << strfrd << endl;
				}
			}
			
			
		}

		if(direction == FR && start + ANCHOR <= refindex[index].ref.length())
		{
			strrev = reverse.substr(kmer_ref[k].second.second, min(ANCHOR, reverse.length() - kmer_ref[k].second.second));
	
			ret_val = find_hband_similarity(refindex[index].ref, strrev, pqueue, start,
				min(refindex[index].ref.length(), end + ANCHOR));
			
			if(ret_val == -1)
			{
				int start1 = kmer_ref[k].second.first;
				int start2 = kmer_ref[k].second.second;
				int kmer_match = 0;
				while(refindex[index].ref.at(start1) == reverse.at(start2))
				{
					kmer_match += 1;
					start1 += 1;
					start2 += 1;
					if(start1 >= refindex[index].ref.length() || start2 >= reverse.length())
						break;
				}
				pre_str_pos = kmer_ref[k].second.second + kmer_match - ANCHOR;
				pre_ref_pos = start + kmer_match - ANCHOR;
				if(pre_str_pos >= 0 && pre_ref_pos >= 1)
				{
					//cout << "max total number of character could be matched = " << kmer_match << endl;
					strrev = reverse.substr(pre_str_pos, min(ANCHOR, reverse.length() - pre_str_pos));
					ret_val = find_hband_similarity(refindex[index].ref, strrev, pqueue, pre_ref_pos,
                                		min(refindex[index].ref.length(), pre_ref_pos + ANCHOR));
					//cout << "backward = " << strrev << endl;
				}	
			}
			
		}

		if(ret_val == -1) {
			end = 0;
			continue;
		}
	
		//cout << "read_pos = " << kmer_ref[k].second.second << ", and ref_pos = " << kmer_ref[k].second.first << endl;
		chain.push_back(make_pair(kmer_ref[k].second.second, kmer_ref[k].second.first));
	}
	
	//t_end = clock();
	//t_anchor += t_end - t_start;
	if(chain.size() < 1) 
	{
		//cout << "The chain size is Zero... " << endl;
		return;
	}


	//cout << "\nIdentifying Max Chain Length\n" << endl;	
	//t_start = clock();

	ref_ind = index;
	read_dir = direction;
	vector<fragment_alignment> optimal_chain;
	int total_len, total_score;

	create_final_alignment(read, read_dir, refindex, ref_ind, chain, optimal_chain, total_len, total_score);  

	//cout << "\nEnd of chain Construction" << endl << endl;
	//t_end = clock();
	//t_extend += t_end - t_start;
	
	if(optimal_chain.size() < 1)
	{
		return;
	}
	

	//t_start = clock();

	int final_gaps, final_mismatches;

	for(int l = 1; l < optimal_chain.size(); l += 2)
	{
		fragment_alignment *fragment_chain;	
		int gaps = 0, mismatches = 0;
		total_len = 0, total_score = 0;
		final_gaps = final_mismatches = 0;
		int current_ref_ind, current_read_ind;
		int total_ref_ind, total_read_ind;
		
		//vector<pair<char, char> > view_chain;
		vector<pair<int, int> > fragment_ind;
		
		total_ref_ind = optimal_chain[l - 1].ref_start;
		total_read_ind = optimal_chain[l - 1].read_start;
		fragment_ind.push_back(make_pair(total_ref_ind, total_read_ind));
			
		fragment_alignment final_alignment_info;

		for(int i = l - 1; i <= l; i++)
		{
			fragment_chain = &optimal_chain[i];
			/*
			cout << "\tfor the index = " << i << ", in direction = " << fragment_chain.read_dir 
					<< ", at reference = " << fragment_chain.ref_ind << endl;
			cout << "\tref_start = " << fragment_chain.ref_start <<
					", and read_start = " << fragment_chain.read_start << endl;
			cout << "\tref_end = " << fragment_chain.ref_end <<
					", and read_end = " << fragment_chain.read_end << endl;
			cout << "\tmatching score = " << fragment_chain.identity_match
					<< " and, alignment len = " << fragment_chain.alignment.size() << ", and ratio = " <<
					(100.000000 * fragment_chain.identity_match / fragment_chain.alignment.size());
			cout << endl;
			*/
			/*
			if(fragment_chain.read_dir == FF)	
				print_alignment(fragment_chain.alignment, refindex[fragment_chain.ref_ind].ref, read,
					fragment_chain.ref_start, fragment_chain.read_start,
					1, fragment_chain, true);
			else
				print_alignment(fragment_chain.alignment, refindex[fragment_chain.ref_ind].ref, reverse,
					fragment_chain.ref_start, fragment_chain.read_start,
					1, fragment_chain, true);
			*/
				
			current_ref_ind = fragment_chain->ref_start;
			current_read_ind = fragment_chain->read_start;
			int index_shifter, index_reverser, newk;
			if(i % 2 == 0)
			{
				index_shifter = fragment_chain->alignment.size() - 1;
				index_reverser = 1;
			}
			else
			{
				index_shifter = 0;
				index_reverser = -1;
			}
	
			for(int k = 0; k < fragment_chain->alignment.size(); k++)//03-20-15
			{
				newk = index_shifter - k * index_reverser;
				if(current_ref_ind == total_ref_ind && current_read_ind == total_read_ind)//added on 01-13-15
				{
				
					if(fragment_chain->alignment[newk].first == fragment_chain->alignment[newk].second
						&& fragment_chain->alignment[newk].first != '-')
					{
						total_score += 1;
					}

					//view_chain.push_back(fragment_chain.alignment[newk]);
					final_alignment_info.alignment.push_back(fragment_chain->alignment[newk]);
					total_len += 1;
					
					if(fragment_chain->alignment[newk].first != '-')
					{
						current_ref_ind += 1;
						total_ref_ind = current_ref_ind;
					}

					if(fragment_chain->alignment[newk].second != '-')
					{
						current_read_ind += 1;
						total_read_ind = current_read_ind;
					}

					if(fragment_chain->alignment[newk].first != fragment_chain->alignment[newk].second
						&& fragment_chain->alignment[newk].first != '-' 
						&& fragment_chain->alignment[newk].second != '-')
							mismatches += 1;

					if(fragment_chain->alignment[newk].first != fragment_chain->alignment[newk].second
						&& (fragment_chain->alignment[newk].first == '-' 
						|| fragment_chain->alignment[newk].second == '-'))
							gaps += 1;

				
				}
				else
				{
					assert(false);//added on 01-13-15
				}
				
				
			}
			/*
			cout << "mark------------------------------------------------------------------------------------" << endl;
			print_alignment(view_chain, result_chain[0].ref_start, result_chain[0].read_start, 1, fragment_chain, true);
			cout << "------------------------------------------------------------------------------------mula" << endl;
			cout << "current_ref_ind = " << current_ref_ind << ", current_read_ind = " << current_read_ind << endl;
			cout << "total_ref_ind = " << total_ref_ind << ", total_read_ind = " << total_read_ind << endl;
			cout << "Total Chain Length = " << total_len << ", and Total Score = " << total_score << endl << endl;
			cout << "mark------------------------------------------------------------------------------------" << endl;
                        print_alignment(view_chain, result_chain[0].ref_start, result_chain[0].read_start, 1, fragment_chain, true);
                        cout << "------------------------------------------------------------------------------------mula" << endl;
			*/
		}

		//cout << "Total Chain Length = " << total_len << ", and Total Score = " << total_score << endl;
		{
			//fragment_alignment final_alignment_info;
			float matching_percentage, coverage_percentage;

			final_alignment_info.total_len = total_len;
			final_alignment_info.identity_match = total_score;
			//final_alignment_info.alignment.clear();
			//final_alignment_info.alignment = view_chain;
			final_alignment_info.ref_start = optimal_chain[l - 1].ref_start;
			final_alignment_info.read_start = optimal_chain[l - 1].read_start;
			final_alignment_info.ref_end = total_ref_ind;//total_ref_ind
			final_alignment_info.read_end = total_read_ind;//total_read_ind
			final_alignment_info.ref_ind = optimal_chain[l - 1].ref_ind;
			final_alignment_info.read_dir = optimal_chain[l - 1].read_dir;
			final_alignment_info.gaps = gaps;
			final_alignment_info.mismatches = mismatches;
			final_alignment_info.fragment_ind = fragment_ind;
			final_gaps = gaps;
			final_mismatches = mismatches;

			if(final_alignment_info.total_len != 0)
				matching_percentage = (100.00 * final_alignment_info.identity_match / final_alignment_info.total_len);
			else
				matching_percentage = 0.0;
                	coverage_percentage = (100.00 * final_alignment_info.total_len / read.length());

			if(matching_percentage > 1.0 * ORIGINALPERCENTMATCH && coverage_percentage > ALIGNMENTCOVERAGE &&
				final_alignment_info.read_end - final_alignment_info.read_start >= MINREADLEN)
			{
				//cout << "RESULT: \t" << (l + 1) << ") match = " << matching_percentage << 
				//	", coverage = " << coverage_percentage << endl;
				list_alignment_info.push_back(final_alignment_info);
			}
			if(DEBUG == 99)
			{	
				
				cout << "ref_start = " << final_alignment_info.ref_start << ", read_start = " << final_alignment_info.read_start
				<< ", ref_end = " << final_alignment_info.ref_end << ", read_end = " << final_alignment_info.read_end << endl;
				cout << "##########################################################################################" << endl << endl;
			
				if(final_alignment_info.read_dir == FF)
					print_alignment(final_alignment_info.alignment, refindex[final_alignment_info.ref_ind].ref, read, 
					final_alignment_info.ref_start, final_alignment_info.read_start, 1, final_alignment_info, true);	
				else
					print_alignment(final_alignment_info.alignment, refindex[final_alignment_info.ref_ind].ref, reverse, 
					final_alignment_info.ref_start, final_alignment_info.read_start, 1, final_alignment_info, true);	
			}
		}
	}
	
	//t_end = clock();
	//t_chain += t_end - t_start;
	
}

void create_final_alignment(string& read, int read_dir, vector<reference_index>& refindex, int ref_ind, vector<pair<int, int> >& chain, 
		vector<fragment_alignment>& fragment_chain, int &total_len, int &total_score) 
{
	vector<pair<int, fragment_alignment> > alignment_list;
	int read_left = 0, ref_left = 0, read_right, ref_right;
	string reverse = reverse_complement(read);
	bool flag_chain[chain.size()], print_flag = false;
	memset(flag_chain, false, sizeof(bool) * chain.size());

	for(int c = 0; c < chain.size(); c++)
	{
		if(flag_chain[c] == true)
			continue;
		/*
		cout << "-----------------------------------------------------------------------------------" << endl;
		cout << "###################################################################################" << endl;
		cout << "-----------------------------------------------------------------------------------" << endl;
		*/
		vector<pair<char, char> > alignment;
		vector<pair<char, char> > current_alignment;
		int read_position = chain[c].first;
		int ref_position = chain[c].second;
		int first_ref_position, last_ref_position;
		int first_read_position, last_read_position;
		int str1_end, str2_end, match;
		string str1, str2, xstring, ystring;

		//cout << "#####extending backward with read_position = " << read_position << 
		//		", and ref_position = " << ref_position << endl;

		read_left = 0;
		ref_left = max(0, chain[c].second - chain[c].first * 1.3);
		read_right = read_position;//+ KMER
		ref_right = ref_position;//+ KMER

		str1 = refindex[ref_ind].ref.substr(ref_left, ref_right - ref_left);	

		if(read_dir == FF)
			str2 = read.substr(read_left, read_right - read_left);
		else
			str2 = reverse.substr(read_left, read_right - read_left);

		reverse_str(str1);
		reverse_str(str2);
		/*
		cout << "str1_len = " << str1.length() << ", and str2_len = " << str2.length() << 
			", min str len = " << min(str1.length(), str2.length()) << endl;
		cout << "string ref  = " << str1.substr(0, min(str1.length(), 60)) << endl;
		cout << "string read = " << str2.substr(0, min(str2.length(), 60)) << endl;
		*/

		fragment_alignment fragment_alignment_backward;

		match = find_banded_similarity(str1, str2, fragment_alignment_backward.alignment, str1_end, str2_end, false);//, 
							//ref_left, read_left, choice, c, -1);//03-26-15
		first_ref_position = ref_right - str1_end;
		first_read_position = read_right - str2_end;
		
		cout << "from banded similarity, str1_end = " << str1_end << ", str2_end = " << str2_end << endl;
		cout << "first_ref_pos = " << first_ref_position << ", and first_read_pos = " << first_read_position << endl;	
		cout << "alignment length = " << alignment.size() << ", with matching score  = " << match << endl;
		
		//for(int a = 0; a < alignment.size(); a++)
   		//	current_alignment.push_back(alignment[a]);
   		//for(int a = alignment.size() - 1; a >= 0;  a--)
		//	current_alignment.push_back(alignment[a]);
		//fragment_alignment fragment_alignment_backward;
		/*
		if(read_dir == FF)
			validate_alignment(refindex, ref_ind, current_alignment, first_ref_position, first_read_position, read);
		else
			validate_alignment(refindex, ref_ind, current_alignment, first_ref_position, first_read_position, reverse);
		*/
		//cout << "Here is the actual alignment = " << endl << endl;
		/*
		if(DEBUG == 99 && print_flag == true)
		{
			if(read_dir == FF)
			{
				print_alignment(current_alignment, refindex[ref_ind].ref, read, first_ref_position, 
				first_read_position, 1, fragment_alignment_backward, false);
			}
			else
			{	
				print_alignment(current_alignment, refindex[ref_ind].ref, reverse, first_ref_position, 
				first_read_position, 1, fragment_alignment_backward, false);
			
			}
			
			assert(fragment_alignment_backward.ref_start == first_ref_position);
			assert(fragment_alignment_backward.ref_end == ref_right - 1);
			assert(fragment_alignment_backward.read_start == first_read_position);
			assert(fragment_alignment_backward.read_end == read_right - 1);
			assert(fragment_alignment_backward.identity_match == match);
			

		}
		else
		*/
		{
			fragment_alignment_backward.ref_start = first_ref_position;
			fragment_alignment_backward.ref_end = ref_right - 1;//<= end
			fragment_alignment_backward.read_start = first_read_position;
			fragment_alignment_backward.read_end = read_right - 1;//<= end
	
			//fragment_alignment_backward.end_to_end = current_alignment;
			//fragment_alignment_backward.alignment = current_alignment;
			fragment_alignment_backward.identity_match = match;
			fragment_alignment_backward.total_len = fragment_alignment_backward.alignment.size();
		}	
		fragment_alignment_backward.read_kmer = read_position;
		fragment_alignment_backward.ref_kmer = ref_position; 
		fragment_alignment_backward.ref_ind = ref_ind;
		fragment_alignment_backward.read_dir = read_dir;
		fragment_chain.push_back(fragment_alignment_backward);
		/*
		cout << "Pushed in the fragment chain with score = " << fragment_alignment_backward.identity_match << 
			", and chain length = " << fragment_chain.size() << endl;
		cout << "\tfragment_read_start = " << fragment_alignment_backward.read_start <<
                        ", fragment_read_end = " << fragment_alignment_backward.read_end << ", read_kmer_position = " <<
                        fragment_alignment_backward.read_kmer << endl;
                cout << "\tfragment_ref_start = " << fragment_alignment_backward.ref_start <<
                        ", fragment_ref_end = " << fragment_alignment_backward.ref_end << ", ref_kmer_position = " <<
                        fragment_alignment_backward.ref_kmer << endl << endl;

		cout << ".................................................................................." << endl;
		cout << "**********************************************************************************" << endl;
		cout << ".................................................................................." << endl;


		cout << "#####extending forward with read_position = " << read_position << 
				", and ref_position = " << ref_position << endl;
		*/
		alignment.clear();
		current_alignment.clear();
		str1_end = str2_end = 0;

		read_left = read_position;//+ KMER
		ref_left = ref_position;// + KMER;
		read_right = read.length();
		ref_right = min(refindex[ref_ind].ref.length(), ref_left + (read.length() - read_left) * 1.3);

		xstring = refindex[ref_ind].ref.substr(ref_left, ref_right - ref_left);  
		if(read_dir == FF)
			ystring = read.substr(read_left, read_right - read_left);
		else
			ystring = reverse.substr(read_left, read_right - read_left);
		/*
		cout << "read_length = " << ystring.length() << ", while ref_length = " << xstring.length() << endl;
		cout << "xstring  = " << xstring.substr(0, min(xstring.length(), 60)) << endl;
		cout << "ystring  = " << ystring.substr(0, min(ystring.length(), 60)) << endl;
		*/

		fragment_alignment fragment_alignment_forward;

		match = find_banded_similarity(xstring, ystring, fragment_alignment_forward.alignment, str1_end, str2_end, true);

		last_ref_position = ref_left + str1_end;
		last_read_position = read_left + str2_end;
		
		cout << "from banded similarity, ref_length = " << str1_end << ", read_length = " << str2_end << endl;
		cout << "last_ref_position = " << last_ref_position << ", and last_read_pos = " << last_read_position << endl;
		cout << "alignment length = " << alignment.size() << ", with matching score  = " << match << endl;
	
		//for(int a = alignment.size() - 1; a >= 0; a--)
		//	current_alignment.push_back(alignment[a]);
		//for(int a = 0; a < alignment.size(); a++)
		//	current_alignment.push_back(alignment[a]);
		//fragment_alignment fragment_alignment_forward;
		/*
		if(read_dir == FF)
			validate_alignment(refindex, ref_ind, current_alignment, ref_left, read_left, read);
		else
			validate_alignment(refindex, ref_ind, current_alignment, ref_left, read_left, reverse);
		*/
		//cout << "Here is the actual alignment = " << endl << endl;
		/*
		if(DEBUG == 99 && print_flag == true)
		{
			if(read_dir == FF)	
			{
				print_alignment(current_alignment, refindex[ref_ind].ref, read, ref_left, 
				read_left, 1, fragment_alignment_forward, false);
			}
			else
			{
				print_alignment(current_alignment, refindex[ref_ind].ref, reverse, ref_left, 
				read_left, 1, fragment_alignment_forward, false);
			}
				
			assert(fragment_alignment_forward.ref_start == ref_left);
			assert(fragment_alignment_forward.ref_end == last_ref_position - 1);
			assert(fragment_alignment_forward.read_start == read_left);
			assert(fragment_alignment_forward.read_end == last_read_position - 1);
			assert(fragment_alignment_forward.identity_match == match);
			
		}
		else
		*/
		{
			fragment_alignment_forward.ref_start = ref_left;
			fragment_alignment_forward.ref_end = last_ref_position - 1;//<= end
			fragment_alignment_forward.read_start = read_left;
			fragment_alignment_forward.read_end = last_read_position - 1;//<= end
	
			//fragment_alignment_forkward.end_to_end = current_alignment;
			//fragment_alignment_forward.alignment = current_alignment;
			fragment_alignment_forward.identity_match = match;
			fragment_alignment_forward.total_len = fragment_alignment_forward.alignment.size();
		}	

		fragment_alignment_forward.read_kmer = read_position;
		fragment_alignment_forward.ref_kmer = ref_position; 
		fragment_alignment_forward.ref_ind = ref_ind;
		fragment_alignment_forward.read_dir = read_dir;
		fragment_chain.push_back(fragment_alignment_forward);
		/*
		cout << "Pushed in the fragment chain with score = " << fragment_alignment_forward.identity_match << 
			", and chain length = " << fragment_chain.size() << endl;
		cout << "\tfragment_read_start = " << fragment_alignment_forward.read_start <<
                        ", fragment_read_end = " << fragment_alignment_forward.read_end << ", read_kmer_position = " <<
                        fragment_alignment_forward.read_kmer << endl;
		cout << "\tfragment_ref_start = " << fragment_alignment_forward.ref_start <<
                        ", fragment_ref_end = " << fragment_alignment_forward.ref_end << ", ref_kmer_position = " <<
                        fragment_alignment_forward.ref_kmer << endl << endl;
		cout << "............................................................................." << endl << endl;
		*/
		for(int x = c + 1; x < chain.size(); x++)//03-26-15
		{
			
			read_position = chain[x].first;
                	ref_position = chain[x].second;

			//cout << "\tcomparing_against: ref_position = " << ref_position << 
			//	", and read_position = " << read_position << endl;

			if(ref_position > fragment_alignment_forward.ref_end)
				break;
			
			if(fragment_alignment_forward.ref_end - ref_position > 0)//FRAGMENTSIZE / 2)
			{
				flag_chain[x] = true;
				//cout << "cancelling the fragment with ref_position = " << ref_position <<
				//	", and read_position = " << read_position << endl;
			}
			
		}
		
	}
	/*
	cout << "Chain Summary: " << endl;
	for(int x = 0; x < fragment_chain.size(); x++)//03-26-15
	{
		cout << "index = " << x << ", fragment_read_start = " << fragment_chain[x].read_start << 
			", fragment_read_end = " << fragment_chain[x].read_end << ", read_kmer_position = " << 
			fragment_chain[x].read_kmer << ", Total length = " << fragment_chain[x].total_len << endl;

		cout << "\tfragment_ref_start = " << fragment_chain[x].ref_start << 
			", fragment_ref_end = " << fragment_chain[x].ref_end << ", ref_kmer_position = " << 
			fragment_chain[x].ref_kmer << ", Total Score = " << fragment_chain[x].identity_match << endl << endl;

	}
	*/
	return;
	/*
	fragment_chain.clear();
	
	cout << "Fragment Chain Summary: " << endl;
	for(int x = 0; x < alignment_list.size(); x++)//03-26-15
	{
		cout << "index = " << alignment_list[x].first << ", fragment_read_start = " << alignment_list[x].second.read_start << 
			", fragment_read_end = " << alignment_list[x].second.read_end << ", read_kmer_position = " << 
			alignment_list[x].second.read_kmer << ", Total length = " << alignment_list[x].second.total_len << endl;

		cout << "\tfragment_ref_start = " << alignment_list[x].second.ref_start << 
			", fragment_ref_end = " << alignment_list[x].second.ref_end << 
			", ref_kmer_position = " << alignment_list[x].second.ref_kmer << 
			", Total Score = " << alignment_list[x].second.identity_match << endl << endl;
		alignment_list[x].second.read_dir = read_dir;
		alignment_list[x].second.ref_ind = ref_ind;
		fragment_chain.push_back(alignment_list[x].second);
	}


	return;
	*/
}


