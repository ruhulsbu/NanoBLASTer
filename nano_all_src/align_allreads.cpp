#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

void align_reads(vector<pair<string, string> >& reference, string& read_file, string& sam_file, vector<reference_index>& refindex)
{
	time_t tstrt, tbgn, tnd;
	time(&tstrt);
		
	/////////////////////////////////////////////////
	ifstream fp_nano;
	string nano_input;
	string nano_file = "last_but_not_nano.txt";
	char *nano = new char[nano_file.length() + 1];
	strcpy(nano, nano_file.c_str());
	fp_nano.open(nano, ifstream::in);
	unordered_map<string, int> nano_read;

	while(getline(fp_nano, nano_input))
	{
		//cout << "" << nano_input << endl;
		nano_read[nano_input] = 1;
		//continue;
	}
	cout << "Total Size of Nano Read = " << nano_read.size() << endl << endl;
	fp_nano.close();
	delete [] nano;
	////////////////////////////////////////////////
	
	ifstream fp_read;
	ofstream fp_sam;
	ofstream fp_2D;
	//vector<string> output;
	
	char *read = new char[read_file.length() + 1];
	strcpy(read, read_file.c_str());

	char *sam = new char[sam_file.length() + 1];
	strcpy(sam, sam_file.c_str());

	//string file_2D_str = logstr.str() + "reverse.fa";
	//string file_2D_str = "2D_file_1K";

	fp_read.open(read, ifstream::in);
	fp_sam.open(sam, ofstream::out | ofstream::app);
	//fp_2D.open(file_2D_str.c_str(), ofstream::out);

	string input, read_name, ref_name;
	string readseq, refgenome;
	string slash = "/";
	int map = 0;
	int count = 0;
	int cant_map = 0;
	int invalid_count = 0;

	fp_csv << "cnt, red_nam, red_len, red_dir, ref_nam, ref_len, ref_pos, score, span, " <<
				"percent, aln_len, spn_rat, aln_tim, tot_tim" << endl;

	while(getline(fp_read, input))
	{
		int find = input.find(slash);
		if(find != string::npos)
			read_name = input.substr(1, find - 1);
		else
			read_name = input.substr(1);

		getline(fp_read, readseq);
		//getline(fp_read, input);
		//getline(fp_read, input);
		//ratio problem with channel_46_read_98_1406145606_2D
		//if(read_name.find("channel_407_read_0_1405093831_2D") == std::string::npos)//to optimize the output
		//if(read_name.find("channel_17_read_24_1405524767_2D") == std::string::npos)//small read to optimize
		//if(read_name.find("channel_201_read_10_1405541481_2D") == std::string::npos)//to compare version 1 and 2
		//if(read_name.find("channel_64_read_7_1403826200_template") == std::string::npos)//max length reads analysis
		//if(read_name.find("channel_424_read_1_1403566249_template") == std::string::npos)//found in last but not in nano
		//if(read_name.find("2D") == std::string::npos)//03-09-2015
		//	continue;//2D== 1D !=
		//if(read_name.find("channel_237_read_42_1406145606_2D") == std::string::npos)
		//if(read_name.find("channel_322_read_11_1405524767_template") == std::string::npos)
		//if(readseq.length() > 100)
		//	continue;
		//if(read_name.find("channel_171_read_2_1403855963_2D") == std::string::npos)//20 times higher than last
		//if(read_name.find("channel_82_read_0_1403855963_2D") == std::string::npos)//20 times higher than last
		//if(read_name.find("channel_221_read_19_1406145606_2D") == std::string::npos)//has maximul length of deletion
		//if(read_name.find("channel_415_read_6_1406242409_template") == std::string::npos)//has 5 times less length than last
		//	continue;
		//if(read_name.find("channel_167_read_19_1403811400_2D") == std::string::npos)//analyze output validity
		//	continue;
		//if(read_name.find("channel_474_read_32_1405524767_template") == std::string::npos)//found in last and nano repeat
		//	continue;
		//if(read_name.find("channel_468_read_12_1403811400_complement") == std::string::npos)//cause exception in nano repeat
		//	continue;
		//if(read_name.find("channel_345_read_7_1403811400_2D") == std::string::npos)//max length increased

		//if(read_name.find("channel_104_read_1_1403551548_template") == std::string::npos)//different in edit not lis
		//if(read_name.find("channel_216_read_0_1403551548_template") == std::string::npos)//different in lis not edit	
		//if(read_name.find("channel_118_read_6_1403551548_template") == std::string::npos)//different in lis and edit	
		//if(read_name.find("channel_486_read_0_1403566249_template") == std::string::npos)//reverse problem
		//if(nano_read.find(read_name) == nano_read.end())
		//	continue;		
		//if(read_name.find("channel_352_read_34_1405541481_template") == std::string::npos)//Why there are multiple results
		//if(read_name.find("channel_68_read_22_1405541481_template") == std::string::npos)//multiple results, boundary problem
		//if(read_name.find("channel_261_read_39_1405541481_template") == std::string::npos)//multiple result indexing
		//	continue;
		//if(read_name.find("channel_302_read_2_1403855963_2D") == std::string::npos)//found in mms not in ssg = align length
		//	continue;
		//if(read_name.find("channel_243_read_0_1403595798_template") == std::string::npos)//found in 40655 not in lis+edit
		//if(read_name.find("channel_452_read_46_1405541481_template") == std::string::npos)//same problem as above
		//	continue;

		//if(read_name.find("channel_431_read_2_1403915857_template") == std::string::npos)//require top 40 tuple list to solve
		//	continue;//solved

		//readseq = readseq.substr(readseq.length() / 2, readseq.length() - readseq.length() / 2);
		//if(read_name.find("channel_199_read_0_1403841073_template") == std::string::npos)//solved
		//if(read_name.find("channel_480_read_91_1406242409_template") == std::string::npos)//in last and not in nano
		//if(read_name.find("channel_389_read_57_1406242409_template") == std::string::npos)//solved
		//	continue;

		//if(read_name.find("channel_56_read_1_1403826200_template") == std::string::npos)
		//	continue;
		//if(read_name.find("channel_356_read_29_1406242409_template") == std::string::npos)// < 40 in nano very weird
		//if(read_name.find("channel_131_read_5_1403826200_template") == string::npos)// < 100 in nano seems weird
		//if(read_name.find("channel_75_read_80_1406145606_template") == std::string::npos)//80% last not found now solved
		//	continue;

		if(count >= MAXREAD && MAXREAD != 0) break;
		count += 1;

		//if(count < 65433) continue;
		cout << count << ") " << read_name << endl;
		if(readseq.length() < MINREADLEN || readseq.length() > MAXREADLEN)//03-09-2015
		{
			cout << "Invalid String Found" << endl;
			invalid_count += 1;
			count -= 1;
			time(&tnd);
			//fp_csv << count << ", " << readseq.length() << ", 0, 0, 0, " <<
			//		"0, 0, 0, 0, 0, 0, 0, " << difftime(tnd, tstrt) << endl;
			continue;
		}
		
		if(count <= MINREAD) continue;
		//if(count < 318) continue;	
	
		time(&tbgn);	
		fp_csv << count << ", " << read_name << ", " << readseq.length() << ", ";

		upper_case(readseq);		
		//reverse_str(readseq);
		//readseq = reverse_complement(readseq);
		
		//fp_2D << input << endl;
		//fp_2D << readseq << endl;
		//continue;

		int match_info, global_match = -1, indpos;
		int match, max_match = 0, match_index, dir;
		vector<vector<string> > list_final_result;
		vector<string> final_result;
		
		time_t start, end;
		clock_t t_start, t_end;

		//for(int i = 0; i < reference.size(); i++)
		{

			//vector<pair<int, pair<int, int> > > kmer_ref;
			vector<pair<int, vector<pair<int, int> > > > kmer_ref;
			cout << "Analysis for forward:" << endl;
			time(&start);
			t_start = clock();

			read_vs_reference(readseq, read_name, FF, refindex, kmer_ref);

			t_end = clock();
			time(&end);
			cout << "Total time taken for calling forward read_vs_ref = " << difftime(end, start) << endl;
			//cout << "Total time for Hash_Lookup = " << (float(t_end - t_start) / CLOCKS_PER_SEC) << endl;
			t_lookup += t_end - t_start;
			//align(readseq, read_name, FF, refindex, kmer_ref, final_result);

			//if(final_result.size() == 0)
			//{
			//kmer_ref.clear();
			cout << "Data for reverse:" << endl;
			time(&start);
			t_start = clock();
			
			string reverse = reverse_complement(readseq);
			//read_vs_reference(reverse, read_name, FR, refindex, kmer_ref);
			
			t_end = clock();
			time(&end);
			cout << "Total time taken for calling reverse read_vs_ref = " << difftime(end, start) << endl;
			//cout << "Total time for Hash_Lookup = " << (float(t_end - t_start) / CLOCKS_PER_SEC) << endl;
			t_lookup += t_end - t_start;
			
			cout << endl <<endl;
			//uncomment here
			list_final_result.clear();
			align(readseq, read_name, FR, refindex, kmer_ref, list_final_result);
			//}
			
			if(list_final_result.size() == 0)
			{
				cant_map += 1;
				kmer_ref.clear();
				time(&tnd);
				//fp_csv << "0, 0, ";
				fp_csv << difftime(tnd, tbgn) << ", " << difftime(tnd, tstrt) << endl;
				//fp_csv << endl;
				continue;
			}
			//assert(final_result.size() == 11);			
			//cout << "No Assertion Occurred for ReferenceIndex = " << read_name << endl;
			//cout << endl << endl;

			kmer_ref.clear();
		}
	
		//cout << "reference_index = " << indpos << ", and match_index = " << match_index <<
		//	", and direction = " << dir << ", with matching = " << max_match << endl;
			
		//break;
		//global_match = read_vs_reference(reference, readseq, read_name, final_result);
		for(int i = 0; i < list_final_result.size(); i++)
		{
			final_result = list_final_result[i];
			fp_sam << final_result[0];
			for(int k = 1; k < final_result.size(); k++)
                	{
                        	fp_sam << "\t" << final_result[k];
                        	//cout << i << ": " << output[k] << endl;
                	}
		
			fp_sam << endl;
			map += 1;
			final_result.clear();
		}
		//if(count >10000)
		//	break;
		
		time(&tnd);
		fp_csv << difftime(tnd, tbgn) << ", " << difftime(tnd, tstrt) << endl;
		cout << "\nTime taken to process " << count << "th read = " << difftime(tnd, tstrt) << "\n" << endl;	
		//break;	
	}

	cout << endl;
	cout << "Overall Statistics - " << endl;
	cout << "total reference size = " << reference.size() << endl << endl;
	cout << "Total read = " << count << endl << endl;
	cout << "Total read mapped = " << map << endl << endl;
	cout << "Total unmapped read = " << cant_map << endl << endl;
	cout << "Out of range read (< 100 or > 15000) = " << invalid_count << endl << endl; 
	//cout << "Total MAX_MATCHED (= " << MAX_MATCHED << ") Reads = " << MAX_SCORED << endl << endl;

	fp_read.close();
	fp_sam.close();
	fp_2D.close();

	delete [] read;
	delete [] sam;

}

bool compare_function(pair<int, int> first_pair, pair<int, int> second_pair)
{
	if(first_pair.first == second_pair.first)
		return (first_pair.second < second_pair.second);
	return (first_pair.first < second_pair.first );
}

void read_vs_reference(string& read, string& read_name, int dir, vector<reference_index>& refindex, 
				vector<pair<int, vector<pair<int, int> > > >& primary_chain)
{
	time_t start, end;

	//unordered_map<long, int> kmer_list;
	//unordered_map<long, int> kmer_map;
	vector<pair<long, int> > kmer_list;
	//vector<pair<int, vector<pair<int, int> > > > primary_chain;
	cout << "\t readseq: " << read_name << " with length = " << read.length() 
			<< " comparing to " << refindex.size() << " references" << endl;
	
	time(&start);

	bool flag = true;
	long hash_key = 0;
	int map_val; 
	int readlen = read.length();
	int prehashval = -1;
	int prehashcnt = 0;

	for(int k = 0; k < read.length(); k++)
	{
		if(flag == true && k + KMER > read.length())
			break;
		for(int l = k, end = k, i = KMER - 2; l < end + KMER - 1 && flag == true; l++, k++, i--)
		{
			map_val = map_value(read.at(l));
			if(map_val == -1)
				break;

			hash_key += base_power_kmer[i][map_val];
			
			//cout << "For character " << read.at(k) << ", at position = " << k <<
			//	", the hash_key = " << hash_key << endl;
		}

		map_val = map_value(read.at(k));
		if(map_val == -1)
		{
			//cout << "Encountered invalid character N ########################" << endl;
			flag = true;
			hash_key = 0;
			continue;
		}
		else
			flag = false;

		hash_key = hash_key * BASE + map_val;
		/*	
		if(kmer_map.find(hash_key) == kmer_map.end())
		{
			kmer_map[hash_key] = 1;//k - KMER + 1;
			kmer_list.push_back(make_pair(hash_key, k - KMER + 1));
		}
		else
		{
			kmer_map[hash_key] += 1;//-1;//need work here
		
			if(kmer_list[kmer_list.size() - 1].first != hash_key)
				kmer_list.push_back(make_pair(hash_key, k - KMER + 1));
			else
				kmer_list[kmer_list.size() - 1].second = k - KMER + 1;
			
		}
		*/
		//cout << "For character " << read.at(k) << ", at position = " << k <<
		//		", the hash_key = " << hash_key << endl;
	 	//if(kmer_map[hash_key] != -1)
	 	
	 	if(prehashval != hash_key)
		{
                	kmer_list.push_back(make_pair(hash_key, k - KMER + 1));
			prehashval = hash_key;
			prehashcnt = 1;
		}
		else
		{
			if(prehashcnt > 2)//suppress same character repeat
                                kmer_list[kmer_list.size() - 1].second = k - KMER + 1;
			else
				kmer_list.push_back(make_pair(hash_key, k - KMER + 1));

			prehashcnt += 1;
		}
		
		map_val = map_value(read.at(k - KMER + 1));
		hash_key -= base_power_kmer[KMER - 1][map_val];
	}

	cout << "Hashing done!!!" << endl;
	int kmer_ref_loc, kmer_read_loc;
	
	for(int index = 0; index < refindex.size(); index++)
	{
		int reflen = refindex[index].ref.length();
		vector<pair<int, int> > kmer_index;

		int interval = 1;
		
		for(int k = 0; k < kmer_list.size(); k++)
		{
			if(refindex[index].index[kmer_list[k].first] == -1)
				continue;
			vector<int> pos = refindex[index].position[refindex[index].index[kmer_list[k].first]];

			for(int i = 0; i < pos.size(); i++)
			{
				kmer_ref_loc = pos[i];
				kmer_read_loc = kmer_list[k].second;

				if(BOUNDED == 0 || kmer_ref_loc - kmer_read_loc / 1.3 >= 0 &&
                                	kmer_ref_loc + (readlen - kmer_read_loc) / 1.3 < reflen)
                                {
                                	kmer_index.push_back(make_pair(kmer_ref_loc, kmer_read_loc));
                                }

			}
		}

		if(kmer_index.size() == 0) 
			continue;

		sort(kmer_index.begin(), kmer_index.end(), compare_function);
		refine_kmer_index(kmer_index, primary_chain, read, dir, refindex, index);
		kmer_index.clear();

	}
	//return;
	
	for(int index = 0; index < refindex.size(); index++)
	{
		int reflen = refindex[index].rev.length();
		vector<pair<int, int> > kmer_index;

		int interval = 1;//introduce constant

		for(int k = 0; k < kmer_list.size(); k++)
		{
			if(refindex[index].revind[kmer_list[k].first] == -1)
				continue;
			vector<int> pos = refindex[index].position[refindex[index].revind[kmer_list[k].first]];

			for(int i = 0; i < pos.size(); i++)
			{
				kmer_ref_loc = reflen - KMER - pos[i];
				kmer_read_loc = readlen - KMER - kmer_list[k].second;

				if(BOUNDED == 0 || kmer_ref_loc - kmer_read_loc / 1.3 >= 0 &&
                                	kmer_ref_loc + (readlen - kmer_read_loc) / 1.3 < reflen)
                                {
                                	kmer_index.push_back(make_pair(kmer_ref_loc, kmer_read_loc));
                                }

			}
		}

		if(kmer_index.size() == 0) 
			continue;

		sort(kmer_index.begin(), kmer_index.end(), compare_function);
		refine_kmer_index(kmer_index, primary_chain, read, FR, refindex, index);
		kmer_index.clear();

	}
	

	time(&end);

	cout << "Total time taken inside read_vs_reference = " << difftime(end, start) << endl;
	cout << "size of the candidate idices = " << primary_chain.size() << endl << endl;

	//kmer_map.clear();	
	kmer_list.clear();
	/*
	sort(primary_chain.begin(), primary_chain.end(), compare_chain);

	for(int i = 0; i < primary_chain.size(); i++)
        {
                cout << "\t" << i << ") " << primary_chain[i].first << " == " << primary_chain[i].second.size() << endl;
                for(int k = 0; k < primary_chain[i].second.size(); k++)
                        cout << "\t\tref = " << primary_chain[i].second[k].first << " : read = " <<
                                primary_chain[i].second[k].second << " : diff = "  <<
                                (primary_chain[i].second[k].first - primary_chain[i].second[k].second) << endl;
        }
	
        primary_chain.clear();
	*/

	return;
}

void create_primary_chain_from_list(node *current_node, vector<pair<int, int> >& chain, node *next_node, node *previous_node,  int readlen)
{
	if(next_node == NULL)
		return;
	//cout << "Chain Ref = " << next_node->ref_ind << ", Read = " << next_node->read_ind << endl;
	//cout << "\t\t\tCurrent Node Ref = " << current_node->ref_ind << ", Read = " << current_node->read_ind << endl;
	bool repeat_sequence = false;
        int range1 = next_node->read_ind - current_node->read_ind;//read
        int range2 = next_node->ref_ind - current_node->ref_ind;//reference

	//cout << "\t\t\trange1 = " << range1 << ", and range2 = " << range2 << endl;
        if(range2 > readlen * 1.30 - current_node->read_ind)// || range2 > 1000)
                return;

        if(range1 < 0)
                create_primary_chain_from_list(current_node, chain, next_node->next, next_node, readlen);
	else if(range2 == 0 || range1 == 0)
	{
		//cout << "\tRepeatDEL: Ref = " << next_node->ref_ind << ", Read = " << next_node->read_ind << endl;
		create_primary_chain_from_list(current_node, chain, next_node->next, next_node, readlen);
		previous_node->next = next_node->next;
		delete next_node;
		next_node = NULL;
	}
        else {
                float diff = (1.0 * range1) / range2;

                if(diff > 0.70 && diff < 1.30)
                {
			//cout << "\tDeleting: Ref = " << next_node->ref_ind << ", Read = " << next_node->read_ind << endl; 
                        if(range1 == range2)
                        {							
				//if(range1 > SEED || repeat_sequence == true)// && chain.size() < 2)//03-09-2015
				{
					chain.push_back(make_pair(next_node->ref_ind, next_node->read_ind));
				}
				
                                //create_primary_chain_from_list(current_node, chain, next_node->next, next_node, readlen);
                                create_primary_chain_from_list(next_node, chain, next_node->next, next_node, readlen);
                                previous_node->next = next_node->next;
				delete next_node;
				next_node = NULL;
                        }
                        else
                        {
                                chain.push_back(make_pair(next_node->ref_ind, next_node->read_ind));

                                create_primary_chain_from_list(next_node, chain, next_node->next, next_node, readlen);
                                previous_node->next = next_node->next;
				delete next_node;
				next_node = NULL;
                        }

                }
                else {
                        create_primary_chain_from_list(current_node, chain, next_node->next, next_node, readlen);
                }
        }

}

void refine_kmer_index(vector<pair<int, int> >& kmer_index, vector<pair<int, vector<pair<int, int> > > >& primary_chain,
                                string read, int dir, vector<reference_index>& refindex, int index)
{
        int first, second, range1, range2;
	int i = 0, k = 0, start, min_len;
        int pre_str_pos, pre_ref_pos;
        int first_index, last_index;
        int ret_val1, ret_val2;
        int min_index = 0;
        int readlen = read.length();
	/*
	vector<pair<int, int> > ref_repeat;
	vector<pair<int, int> > read_repeat;

        for(i = 0; i < kmer_index.size(); i++)
        {
                if(kmer_index[i].first == -1)
                        continue;

                vector<pair<int, int> > chain;
                chain.push_back(kmer_index[i]);

                create_primary_chain(kmer_index, chain, i, i + 1, readlen);
		if(chain.size() > 1)
                	primary_chain.push_back(make_pair(index * MAXLEN + dir, chain));

        }
	return;
	*/

	node *head = new node;
	head->ref_ind = kmer_index[i].first;
	head->read_ind = kmer_index[i].second;
	head->next = NULL;
	node *current_node = head;

	//cout << "\nIndex = " << index << " ********************Starting New Analysis******************** = " 
	//		<< kmer_index.size() << endl << endl;

	int ref_dif, read_dif;
	int next_ref_dif, next_read_dif;
	int prev_ref_dif, prev_read_dif;

	bool ref_repeat_flag = false;
	bool read_repeat_flag = false;
	int best_index, max_length = 0;
	int ref_start, read_start, l, x;
	
	for(i = 1; i < kmer_index.size(); i++)
	{
				
 		//cout << i << ") Current Ref = " << current_node->ref_ind << ", Read = " << current_node->read_ind << endl;
		ref_dif = (kmer_index[i].first - current_node->ref_ind);
		read_dif = (kmer_index[i].second - current_node->read_ind);

		node *next_node = new node;
		next_node->ref_ind = kmer_index[i].first;
		next_node->read_ind = kmer_index[i].second;
		next_node->next = NULL;
		//cout << "Considering " << k << ") Reference = " << next_node->ref_ind << ", Read = " << next_node->read_ind << endl;
		current_node->next = next_node;
		current_node = next_node;
		k += 1;
		max_length = 0;
		ref_repeat_flag = false;
		read_repeat_flag = false;
	}

	if(i < kmer_index.size())
	{
		node *next_node = new node;
		next_node->ref_ind = kmer_index[i].first;
		next_node->read_ind = kmer_index[i].second;
		next_node->next = NULL;
		//cout << "\t" << k << ") Reference = " << next_node->ref_ind << ", Read = " << next_node->read_ind << endl;
		current_node->next = next_node;
		current_node = next_node;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	node *tmp_node;
	current_node = head;
	
	while(current_node != NULL)
	{
		//cout << "\nStarting Chain Creation: " << endl;
		//cout << "Ref = " << current_node->ref_ind << ", And Read = " << current_node->read_ind << endl;
		//if(current_node->read_ind < readlen / 2)
		{
			vector<pair<int, int> > chain;
                	chain.push_back(make_pair(current_node->ref_ind, current_node->read_ind));

                	create_primary_chain_from_list(current_node, chain, current_node->next, current_node, readlen);
			if(chain.size() > 1)
			{
                		primary_chain.push_back(make_pair(index * MAXLEN + dir, chain));
				//cout << "Size of Chain is = " << chain.size() << endl; 
			}
		}
		
		tmp_node = current_node;
		current_node = current_node->next;
		delete tmp_node;
		tmp_node = NULL;
	
		//cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl << endl;
	}


	/*		
	current_node = head;
	while(current_node != NULL)
	{
		tmp_node = current_node;
		cout << "Exception Reference = " << tmp_node->ref_ind << ", Read = " << tmp_node->read_ind << endl;
		current_node = current_node->next;
		delete tmp_node;
		tmp_node = NULL;
	}
	head = NULL;

	
	tmp_node = head;
	while(tmp_node != NULL)
	{
		cout << "Extra Reference = " << tmp_node->ref_ind << ", Read = " << tmp_node->read_ind << endl;
		tmp_node = tmp_node->next;
	}
	*/
}

void create_primary_chain(vector<pair<int, int> >& kmer_index, vector<pair<int, int> >& chain,
                int start_ind, int next_ind, int readlen)
{
	if(next_ind >= kmer_index.size())
		return;

        int range1 = kmer_index[next_ind].second - kmer_index[start_ind].second;//read
        int range2 = kmer_index[next_ind].first - kmer_index[start_ind].first;//reference

        if(range2 > readlen * 1.30 - kmer_index[start_ind].second)
                return;

        if(kmer_index[next_ind].first == -1 || range1 < 0)
                create_primary_chain(kmer_index, chain, start_ind, next_ind + 1, readlen);
        else {
                float diff = (1.0 * range1) / range2;

                if(diff > 0.70 && diff < 1.30)
                {
                        if(range1 == range2)
                        {
							
				if(range1 > SEED && chain.size() < 2)//03-09-2015
				{
					chain.push_back(kmer_index[next_ind]);
					//start_ind = next_ind;
				}
				
                                create_primary_chain(kmer_index, chain, start_ind, next_ind + 1, readlen);
                                kmer_index[next_ind] = make_pair(-1, -1);
                        }
                        else
                        {
                                chain.push_back(kmer_index[next_ind]);
                                start_ind = next_ind;

                                create_primary_chain(kmer_index, chain, start_ind, next_ind + 1, readlen);
                                kmer_index[next_ind] = make_pair(-1, -1);
                        }
                }
                else {
                        create_primary_chain(kmer_index, chain, start_ind, next_ind + 1, readlen);
                }
        }
}


