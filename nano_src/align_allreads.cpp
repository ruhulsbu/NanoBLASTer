#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"


void align_reads(vector<pair<string, string> >& reference, string& read_file, string& sam_file, vector<reference_index>& refindex)
{
	time_t tstrt, tbgn, tnd;
	time(&tstrt);
		
	/*//Can be used for analyzing the difference between LAST and NanoBLASTer
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
	*/
	
	ifstream fp_read;
	ofstream fp_sam;
	
	char *read = new char[read_file.length() + 1];
	strcpy(read, read_file.c_str());

	char *sam = new char[sam_file.length() + 1];
	strcpy(sam, sam_file.c_str());

	fp_read.open(read, ifstream::in);
	fp_sam.open(sam, ofstream::out | ofstream::app);

	string input, read_name, ref_name;
	string readseq, refgenome;
	string slash = "/";
	int map = 0;
	int count = 0;
	int cant_map = 0;
	int invalid_count = 0;

	fp_csv << "cnt, red_nam, red_len, red_dir, ref_nam, ref_len, ref_pos, score, span, " <<
				"percent, aln_len, spn_rat, aln_tim, tot_tim" << endl;

	getline(fp_read, input);
	while(!fp_read.eof())
	{
		//int find = input.find(slash);
		//if(find != string::npos)
		//	read_name = input.substr(1, find - 1);
		//else
		//	read_name = input.substr(1);
		read_name = input.substr(1);

		readseq = "";
		while(getline(fp_read, input))
		{
			if(input.length() == 0)
				continue;
			if(input.at(0) == '>')
				break;

			readseq += input;
		}
		
		//getline(fp_read, input);
		//getline(fp_read, input);

		//ratio problem with channel_46_read_98_1406145606_2D
		//if(read_name.find("channel_407_read_0_1405093831_2D") == std::string::npos)//to optimize the output
		//if(read_name.find("channel_17_read_24_1405524767_2D") == std::string::npos)//small read to optimize
		//if(read_name.find("channel_201_read_10_1405541481_2D") == std::string::npos)//to compare version 1 and 2
		//if(read_name.find("channel_64_read_7_1403826200_template") == std::string::npos)//max length reads analysis
		//if(read_name.find("channel_424_read_1_1403566249_template") == std::string::npos)//found in last but not in nano
		//if(read_name.find("2D") == std::string::npos)//03-09-2015
		//if(read_name.find("channel_237_read_42_1406145606_2D") == std::string::npos)
		//if(read_name.find("channel_322_read_11_1405524767_template") == std::string::npos)

		//if(read_name.find("channel_171_read_2_1403855963_2D") == std::string::npos)//20 times higher than last
		//if(read_name.find("channel_82_read_0_1403855963_2D") == std::string::npos)//20 times higher than last
		//if(read_name.find("channel_221_read_19_1406145606_2D") == std::string::npos)//has maximul length of deletion
		//if(read_name.find("channel_415_read_6_1406242409_template") == std::string::npos)//has 5 times less length than last
		//if(read_name.find("channel_167_read_19_1403811400_2D") == std::string::npos)//analyze output validity
		//if(read_name.find("channel_474_read_32_1405524767_template") == std::string::npos)//found in last and nano repeat
		//if(read_name.find("channel_468_read_12_1403811400_complement") == std::string::npos)//cause exception in nano repeat
		//if(read_name.find("channel_345_read_7_1403811400_2D") == std::string::npos)//max length increased

		//if(read_name.find("channel_104_read_1_1403551548_template") == std::string::npos)//different in edit not lis
		//if(read_name.find("channel_216_read_0_1403551548_template") == std::string::npos)//different in lis not edit	
		//if(read_name.find("channel_118_read_6_1403551548_template") == std::string::npos)//different in lis and edit	
		//if(read_name.find("channel_486_read_0_1403566249_template") == std::string::npos)//reverse problem
		//if(nano_read.find(read_name) == nano_read.end())
		//if(read_name.find("channel_352_read_34_1405541481_template") == std::string::npos)//Why there are multiple results
		//if(read_name.find("channel_68_read_22_1405541481_template") == std::string::npos)//multiple results, boundary problem
		//if(read_name.find("channel_261_read_39_1405541481_template") == std::string::npos)//multiple result indexing

		//if(read_name.find("channel_302_read_2_1403855963_2D") == std::string::npos)//found in mms not in ssg = align length
		//if(read_name.find("channel_243_read_0_1403595798_template") == std::string::npos)//found in 40655 not in lis+edit
		//if(read_name.find("channel_452_read_46_1405541481_template") == std::string::npos)//same problem as above
		//if(read_name.find("channel_431_read_2_1403915857_template") == std::string::npos)//require top 40 tuple list to solve
		//readseq = readseq.substr(readseq.length() / 2, readseq.length() - readseq.length() / 2);
		//if(read_name.find("channel_199_read_0_1403841073_template") == std::string::npos)//solved
		//if(read_name.find("channel_480_read_91_1406242409_template") == std::string::npos)//in last and not in nano
		//if(read_name.find("channel_389_read_57_1406242409_template") == std::string::npos)//solved
		//if(read_name.find("channel_56_read_1_1403826200_template") == std::string::npos)
		//if(read_name.find("channel_356_read_29_1406242409_template") == std::string::npos)// < 40 in nano very weird
		//if(read_name.find("channel_131_read_5_1403826200_template") == string::npos)// < 100 in nano seems weird
		//if(read_name.find("channel_75_read_80_1406145606_template") == std::string::npos)//80% last not found now solved
		//	continue;


		if(count >= MAXREAD && MAXREAD != 0) break;
		count += 1;

		//if(count < 11762) continue;
		cout << count << ") " << read_name << endl;
		if(readseq.length() < MINREADLEN || readseq.length() > MAXREADLEN)//03-09-2015
		{
			cout << "Invalid String Found" << endl;
			invalid_count += 1;
			count -= 1;
			fp_sam << read_name << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << readseq << "\t*" << endl;
			time(&tnd);
			//fp_csv << count << ", " << readseq.length() << ", 0, 0, 0, " <<
			//		"0, 0, 0, 0, 0, 0, 0, " << difftime(tnd, tstrt) << endl;
			continue;
		}
		
		if(count <= MINREAD) continue;
		//if(count < 318) continue;	
	
		time(&tbgn);
		if(DEBUG == 99)	
			fp_csv << count << ", " << read_name << ", " << readseq.length() << ", ";

		upper_case(readseq);		
		//reverse_str(readseq);
		//readseq = reverse_complement(readseq);
		
		int match_info, global_match = -1, indpos;
		int match, max_match = 0, match_index, dir;
		vector<vector<string> > list_final_result;
		
		//time_t start, end;
		//clock_t t_start, t_end;

		//for(int i = 0; i < reference.size(); i++)
		{

			//vector<pair<int, pair<int, int> > > kmer_ref;
			vector<pair<int, vector<pair<int, int> > > > kmer_ref;
			//cout << "Analysis for forward:" << endl;
			//
			//time(&start);
			//t_start = clock();

			read_vs_reference(readseq, read_name, FF, refindex, kmer_ref);

			//t_end = clock();
			//time(&end);
			//cout << "Total time taken for calling forward read_vs_ref = " << difftime(end, start) << endl;
			
			//t_lookup += t_end - t_start;
			//align(readseq, read_name, FF, refindex, kmer_ref, final_result);

			//cout << "Data for reverse:" << endl;
			//time(&start);
			//t_start = clock();
		
			if(SINGLE == 1)
			{	
				string reverse = reverse_complement(readseq);
				read_vs_reference(reverse, read_name, FR, refindex, kmer_ref);
			}
			//t_end = clock();
			//time(&end);
			//cout << "Total time taken for calling reverse read_vs_ref = " << difftime(end, start) << endl;
			
			//t_lookup += t_end - t_start;
			//cout << endl <<endl;

			//uncomment here for aligninng read
			list_final_result.clear();
			align(readseq, read_name, FR, refindex, kmer_ref, list_final_result);
			
			
			if(list_final_result.size() == 0)
			{
				cant_map += 1;
				kmer_ref.clear();
				//time(&tnd);
				//fp_csv << difftime(tnd, tbgn) << ", " << difftime(tnd, tstrt) << endl;
				fp_sam << read_name << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << readseq << "\t*" << endl;
				continue;
			}

			kmer_ref.clear();
		}
	
		for(int i = 0; i < list_final_result.size(); i++)
		{
			vector<string>& final_result = list_final_result[i];
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
		/*
		if(list_final_result.size() == 0 && SAM_FORMAT == 1)
		{
			fp_sam << read_name << "\t4\t*\t0\t0\t*\t*\t0\t0\t" << readseq << "\t*" << endl;
		}
		*/
		//time(&tnd);
		list_final_result.clear();
		if(DEBUG == 99)
			fp_csv << endl;//difftime(tnd, tbgn) << ", " << difftime(tnd, tstrt) << endl;
		//cout << "\nTime taken to process " << count << "th read = " << difftime(tnd, tstrt) << "\n" << endl;	
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
	vector<pair<long, int> > kmer_list;
	//unordered_map<long, int> kmer_list;
	//unordered_map<long, int> kmer_map;
	//vector<pair<int, vector<pair<int, int> > > > primary_chain;

	//cout << "\t readseq: " << read_name << " with length = " << read.length() 
	//		<< " comparing to " << refindex.size() << " references" << endl;
	
	time(&start);

	bool flag = true;
	long hash_key = 0;
	int map_val; 
	int readlen = read.length();
	long prehashval = -1;
	int prehashcnt = 0;

	for(int k = 0; k < read.length(); k++)
	{
		if(flag == true)
		{
			if(k + KMER > read.length())
				break;
			for(int l = k, end = k, i = KMER - 2; l < end + KMER - 1; l++, k++, i--)
			{
				map_val = map_value(read.at(l));
				if(map_val == -1)
					break;

				hash_key += base_power_kmer[i][map_val];
			
				//cout << "For character " << read.at(k) << ", at position = " << k <<
				//	", the hash_key = " << hash_key << endl;
			}
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
		//cout << "For character " << read.at(k) << ", at position = " << k <<
		//		", the hash_key = " << hash_key << endl;
	 	*/
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

	{
		cout << "Starting KMER Chain Analysis for " << dir << endl;
		vector<vector<pair<int, int> > > kmer_index;
		for(int i = 0; i < refindex.size(); i++)
		{
			vector<pair<int, int> > kmer_pair;
			kmer_index.push_back(kmer_pair);
		}

		int interval = 1, index, reflen;
		int kmer_ref_loc, kmer_read_loc;

		for(int k = 0; k < kmer_list.size(); k++)
		{
			if(basic_index.index[kmer_list[k].first] == -1)
				continue;
			//vector<int> pos = basic_index.position[basic_index.index[kmer_list[k].first]];
			int ref_location;
			//for(int i = 0; i < pos.size(); i++)
			for(vector<int>::iterator it = basic_index.position[basic_index.index[kmer_list[k].first]].begin(); 
				it != basic_index.position[basic_index.index[kmer_list[k].first]].end(); ++it) 
			{
				ref_location = *it;
				if(ref_location < 0)//if(pos[i] < 0)
				{
					index = abs(ref_location) - 1;//abs(pos[i]) - 1;
					reflen = refindex[index].ref.length();
					continue;
				}
				kmer_ref_loc = ref_location;// pos[i];
				kmer_read_loc = kmer_list[k].second;
				if(CIRCULAR == 0 || kmer_ref_loc - kmer_read_loc / 1.3 >= 0 &&
					kmer_ref_loc + (readlen - kmer_read_loc) / 1.3 < reflen)
				{
					//cout << "index = " << index << endl;
					//assert(index >= 0 && index < refindex.size());
					kmer_index[index].push_back(make_pair(kmer_ref_loc, kmer_read_loc));
				}
			}
		}

		for(int index = 0; index < refindex.size(); index++)
		{
			if(kmer_index[index].size() == 0)
				continue;
			sort(kmer_index[index].begin(), kmer_index[index].end(), compare_function);
			/*for(int i = 0; i < kmer_index[index].size(); i++)
			{
				cout << "first = " << kmer_index[index][i].first << ", and second = " << kmer_index[index][i].second << endl;
			}*/
			refine_kmer_index(kmer_index[index], primary_chain, read, dir, refindex, index);
			kmer_index[index].clear();
		}
	}

	if(SINGLE == 1)
	{
		kmer_list.clear();
		return;
	}

	{
		cout << "Starting Reverse Analysis for " << FR << endl;
		vector<vector<pair<int, int> > > kmer_index;
		kmer_index.clear();
		for(int i = 0; i < refindex.size(); i++)
		{
			vector<pair<int, int> > kmer_pair;
			kmer_index.push_back(kmer_pair);
		}

		int interval = 1, index, reflen;
		int kmer_ref_loc, kmer_read_loc;

		for(int k = 0; k < kmer_list.size(); k++)
		{
			if(basic_index.revind[kmer_list[k].first] == -1)
				continue;
			//vector<int> pos = basic_index.position[basic_index.revind[kmer_list[k].first]];
			int ref_location;
			//for(int i = 0; i < pos.size(); i++)
			for(vector<int>::iterator it = basic_index.position[basic_index.revind[kmer_list[k].first]].begin();
                                it != basic_index.position[basic_index.revind[kmer_list[k].first]].end(); ++it)
			{
				ref_location = *it;
				if(ref_location < 0)//if(pos[i] < 0)
				{
					index = abs(ref_location) - 1;//index = abs(pos[i]) - 1;
					reflen = refindex[index].rev.length();
					continue;
				}
				kmer_ref_loc = reflen - KMER - ref_location;
				//kmer_ref_loc = reflen - KMER - pos[i];
				kmer_read_loc = readlen - KMER - kmer_list[k].second;
				if(CIRCULAR == 0 || kmer_ref_loc - kmer_read_loc / 1.3 >= 0 &&
					kmer_ref_loc + (readlen - kmer_read_loc) / 1.3 < reflen)
				{
					//cout << "index = " << index << endl;
					//assert(index >= 0 && index < refindex.size());
					kmer_index[index].push_back(make_pair(kmer_ref_loc, kmer_read_loc));
				}
			}
		}

		for(int index = 0; index < refindex.size(); index++)
		{
			if(kmer_index[index].size() == 0)
				continue;
			sort(kmer_index[index].begin(), kmer_index[index].end(), compare_function);
			/*for(int i = 0; i < kmer_index[index].size(); i++)
			{
				cout << "first = " << kmer_index[index][i].first << ", and second = " << kmer_index[index][i].second << endl;
			}*/
			refine_kmer_index(kmer_index[index], primary_chain, read, FR, refindex, index);
			kmer_index[index].clear();
		}
	}



	time(&end);
	//cout << "Total time taken inside read_vs_reference = " << difftime(end, start) << endl;
	//cout << "size of the candidate idices = " << primary_chain.size() << endl << endl;

	//kmer_map.clear();	
	kmer_list.clear();
	
	return;
}

void create_primary_chain_from_list(node *head_node, vector<pair<int, int> >& chain, node *next_node, node *previous_node,  int readlen)
{
        bool repeat_sequence;
        int range1, range2;
	int max_reference_delta;
	
	node *current_node = new node;
        current_node->ref_ind = head_node->ref_ind;
        current_node->read_ind = head_node->read_ind;
        current_node->next = head_node->next;

        while(next_node != NULL)
        {
                //cout << "Chain Ref = " << next_node->ref_ind << ", Read = " << next_node->read_ind << endl;
                //cout << "\t\t\tCurrent Node Ref = " << current_node->ref_ind << ", Read = " << current_node->read_ind << endl;
                repeat_sequence = false;
                range1 = next_node->read_ind - current_node->read_ind;//read
                range2 = next_node->ref_ind - current_node->ref_ind;//reference

		max_reference_delta = min(1.30 * readlen - current_node->read_ind, SEEDTUPLEDIST);
		if(range2 > max_reference_delta)
                //if(range2 > readlen * 1.30 - current_node->read_ind)// || range2 > 1000)
		{
                        break;
		}
                if(range1 < 0)
                {
                        previous_node = next_node;
                        next_node = next_node->next;
                        continue;
                }
                else if(range2 == 0 || range1 == 0)
                {
                        //create_primary_chain_from_list(current_node, chain, next_node->next, next_node, readlen);
                        previous_node->next = next_node->next;
                        delete next_node;
                        next_node = previous_node->next;//NULL;
                        //return;
                }
                else {
                        float diff = (1.0 * range1) / range2;

			//if(diff > 0.90 && diff < 1.10) // hack
                        if(diff > 0.70 && diff < 1.30)
                        {
				current_node->ref_ind = next_node->ref_ind;
                                current_node->read_ind = next_node->read_ind;
                                current_node->next = next_node->next;
				
                                //if(range1 == range2)
                                {
					chain.push_back(make_pair(next_node->ref_ind, next_node->read_ind));
					//cout << "Chain Ref = " << next_node->ref_ind << ", Read = " << next_node->read_ind << endl;
                                        //create_primary_chain_from_list(next_node, chain, next_node->next, next_node, readlen);
                                        previous_node->next = next_node->next;
                                        delete next_node;
                                        next_node = NULL;
                                }
                                /*else
                                {
                                        chain.push_back(make_pair(next_node->ref_ind, next_node->read_ind));

                                        //create_primary_chain_from_list(next_node, chain, next_node->next, next_node, readlen);
                                        previous_node->next = next_node->next;
                                        delete next_node;
                                        next_node = NULL;
                                }
				*/
				next_node = previous_node->next;
                                //return;
                        }
                        else {
                                previous_node = next_node;
                                next_node = next_node->next;
                        }
                }
        }

	delete current_node;
	current_node = NULL;
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
                		primary_chain.push_back(make_pair(index * MULTIPLIER + dir, chain));
				//cout << "Size of Chain is = " << chain.size() << endl; 
			}
		}
		
		tmp_node = current_node;
		current_node = current_node->next;
		delete tmp_node;
		tmp_node = NULL;
	
		//cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl << endl;
	}

}


