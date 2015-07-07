
void find_seed(string& str1, string& str2, priority_queue<long>& pqueue)
{
	unordered_map<string, int> map;
	priority_queue<long> temporary;
	vector<long> enlist;

	int i, k, length, hband, index;
	int match, max_local_match;
	int first = 0, second = 0, count = 0;
	long placement, preindex, postindex;
	double score;
	string kmer;

	hband = min(2 * HBAND, str2.length());
	if(DEBUG == 6)
		cout << "Size of HBAND = " << hband << " and WINDOW = " << WINDOW << endl;

	for(i = 0; i < hband - WINDOW + 1; i++)
	{
		kmer = str2.substr(i, WINDOW);
		if(map.find(kmer) == map.end())
		{
			map[kmer] = 1;
			first += 1;
		}
		/*else 
			map[kmer] += 1;
		*/
	}

	preindex = 0;
	for(i = 0; i < hband - WINDOW + 1; i++)
	{
		kmer = str1.substr(i, WINDOW);
		if(map.find(kmer) != map.end())
		{
			if(map[kmer] == 1)
			{
				map[kmer] -= 1;
				second += 1;
				if(DEBUG == 7) cout << "KMER MATCH = " << kmer << endl;
				/*	
				if(preindex + WINDOW > i)
					count += i - preindex;
				else
					count += WINDOW;

				preindex = i;
				*/
			}
		}
	}

	if(DEBUG == 7) cout << "first = " << first << " and second = " << second << endl;
	index = 1;
	score = second / first;
	placement = second * MAXLEN + (MAXLEN - index);
	assert(placement >= 0);
	enlist.push_back(placement);
	//temporary.push(placement);
	if(DEBUG == 7) cout << "index = " << index << " and matching = " << second 
					<< " with placement = " << placement << endl;

	//if(DEBUG == 6) cout << "starting placement = " << placement << endl;
	/*
	if(second > 0)
		cout << "Matching at index = 0 is " << second << endl; 
	*/
	max_local_match = 0;
	for(k = 0; i < str1.length() - WINDOW + 1; k++, i++)
	{
		kmer = str1.substr(k, WINDOW);
		if(map.find(kmer) != map.end())
		{
			if(map[kmer] == 0)
			{
				map[kmer] += 1;
				second -= 1;
				if(DEBUG == 7) cout << "DOING KMER DELETION = " << kmer << endl;
			}
		}

		kmer = str1.substr(i, WINDOW);
		if(map.find(kmer) != map.end())
		{
			if(map[kmer] == 1)
			{
				map[kmer] -= 1;
				second += 1;
				if(DEBUG == 7) cout << "KMER MATCH = " << kmer << endl;
			}
		}
		
		index = i + 1 - (hband - WINDOW);
		score = second / first;
		placement = second * MAXLEN + (MAXLEN - index);
		assert(placement >= 0);
		if(DEBUG == 7) cout << "index = " << index << " and matching = " << second 
					<< " with placement = " << placement << endl;
		//temporary.push(placement); //08.19.2014
		if(max_local_match < second)
			max_local_match = second;
		enlist.push_back(placement);
		/*
		if(second > 0)
			cout << "Matching at index = " << k << " is " << second << endl;
		*/
	}
	
	for(k = 0; k < enlist.size(); k++)
	{
		index = MAXLEN - enlist[k] % MAXLEN;
		match = enlist[k] / MAXLEN;

		if(match < floor(max_local_match * HRATIO))
			continue;
		placement = max_local_match * MAXLEN + (MAXLEN - index);

		if(DEBUG == 7) cout << "index = " << index << " and matching = " << match 
					<< " with placement = " << placement << endl;
		temporary.push(placement);
	}

	if(max_local_match < hband / 200 + 1)//CHANGE
		return;
	
	if(DEBUG == 7)
		cout << "Size of the Temporary Queue = " << temporary.size() << "" << endl; 
	//pqueue = temporary;
	//return;

	preindex = MAXLEN - (temporary.top() % MAXLEN);
	max_local_match = temporary.top() / MAXLEN;
	postindex = preindex;
	bool flag = true;

	while(!temporary.empty())
	{
		index = MAXLEN - temporary.top() % MAXLEN;
		match = temporary.top() / MAXLEN;

		if(DEBUG == 7) cout << "---------------------------------Index = " << 
			index << " with Match = " << match << " and Value = " << temporary.top() << endl;
		
		//if(index >= preindex && index - preindex <= hband
		if(index >= postindex && index - postindex <= hband
			&& match >= floor(max_local_match * HRATIO))
		{	
			postindex = index;
			temporary.pop();
			continue;
		}

		if(preindex >= index && preindex - index <= hband &&
			match >= floor(max_local_match * HRATIO))
		{
			preindex = index;
			temporary.pop();
			continue;
		}

		if(index >= preindex && index <= postindex)
		{
			temporary.pop();
			continue;
		}
		
		if(DEBUG == 6)
			cout << "Calling Hband with " << preindex << " And " << (postindex + hband) << endl;
		time_t tstart, tend;
		time(&tstart);
		find_hband_similarity(str1, str2, pqueue, preindex, 
					min(postindex + hband, str1.length()));	
		time(&tend);
		cout << "==================== Timke taken in find_seed for find_hband_sim = " << difftime(tend, tstart) << endl;
	
		preindex = postindex = index;

		temporary.pop();
		if(match < floor(max_local_match * HRATIO))
		{
			flag = false;
			break;
		}

	}

	if(pqueue.size() == 0 || flag == true)
	{
		if(DEBUG == 6)
			cout << "Calling Hband with " << preindex << " And " << min(postindex + hband, str1.length()) << endl;

		find_hband_similarity(str1, str2, pqueue, preindex, min(postindex + hband, str1.length()));
	}

	if(DEBUG == 6) cout << "PQUEUE Size = " << pqueue.size() << endl;	
	map.clear();
}

bool compare(const long first, const long second)
{
	long index1 = MAXLEN - first % MAXLEN;
	long index2 = MAXLEN - second % MAXLEN;

	return (index1 < index2);
}

long neighbor_search(vector<long>& indexvec, long lookfor)
{
	long placement, match, index;
	int begin = 0, middle, end = indexvec.size() - 1;
	//int middle = (begin + end) / 2;
	long difference = MAXLEN, closest;
	//cout << "Searching = " << lookfor << endl; 
	while(begin <= end)
	{	
		middle = (begin + end) / 2;
		placement = indexvec[middle];
		index = MAXLEN - placement % MAXLEN;
		//cout << "Hit At = " << index << endl;
		if(index == lookfor)
			break;

		if(lookfor < index)
			end = middle - 1;
		else
			begin = middle + 1;
		
	}

	return indexvec[middle];
}

void randomized_hband(string& str1, string& str2, priority_queue<long>& pqueue)
{
	priority_queue<long> xqueue;
	priority_queue<long> yqueue;
	vector<long> indexvec;

	time_t tstart, tend;
	int hband = min(2 * HBAND, str2.length() / 2);
	int count = 0;

	string begin = str2.substr(0, hband);
	string end = str2.substr(str2.length() - hband);
	time(&tstart);
	find_seed(str1, begin, xqueue);
	time(&tend);
	cout << "==================== Timke taken in randomized_hband for find_seed = " << difftime(tend, tstart) << endl;
	time(&tstart);
	find_seed(str1, end, yqueue);
	time(&tend);
	cout << "==================== Timke taken in randomized_hband for find_seed = " << difftime(tend, tstart) << endl;
	
	long first, second, lookfor, placement;
	long approximate_len;
	while(!yqueue.empty())
	{
		indexvec.push_back(yqueue.top());
		yqueue.pop();
	}
	sort(indexvec.begin(), indexvec.end(), compare);
	/*	
	for(vector<long>::iterator it = indexvec.begin(); it != indexvec.end(); it++)
	{
		cout << "sorted element for index of end part: " << *it << endl;
	}
	*/

	while(!xqueue.empty() && count < MAXTRY)
	{
		first = MAXLEN - xqueue.top() % MAXLEN;
		lookfor = first + str2.length() - hband;
		//cout << "Chudir Vair Point !!!" << endl;		
		placement = neighbor_search(indexvec, lookfor);
		second = MAXLEN - placement % MAXLEN;
		approximate_len = second - first + hband;
		if(approximate_len < (str2.length() * (1 + ERROR)) && approximate_len > (str2.length() * (1 - ERROR)))
		{
			pqueue.push(xqueue.top()); 
			count += 1;
			cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
		}
		cout << "Begining index = " << first << " and Searched for index " << lookfor << " And got = " << second << endl;
		xqueue.pop();
	}
	/*	
	cout << "WHATS UP!"  << endl;
	
	while(!xqueue.empty() && !yqueue.empty())
	{
		first = MAXLEN - xqueue.top() % MAXLEN;
		second = MAXLEN - yqueue.top() % MAXLEN;
		cout << "First = " << first << ", AND Second = " << second << endl;

		xqueue.pop();
		yqueue.pop();
	}
	*/
}

int is_overlap(string& str1, string& str2, SeqPQ& sequences, int dir)
{
	int i, row, column;
	int index, match, count = 0;
	int max_local_match;
	long dequeue_value;//, **path = NULL;
	float kband_ratio;
	int local_index, global_match = -MAXLEN;
	string xstring, ystring;
	priority_queue<long> pqueue;
	vector<pair<char, char> > alignment, max_alignment;
	
	time_t start, end;

	/*	
	local_index = find_hband_similarity(str1, str2, pqueue);

	if(DEBUG == 5)	printf("The line number is (completed hband): %d\n", __LINE__);

	if(local_index < 0)
		return -1;
	*/
	//cout << str1 << endl;
	//cout << str2 << " In Direction = " << dir << endl;
	
	time(&start);
	randomized_hband(str1, str2, pqueue);
	time(&end);
	cout << "==================== Timke taken in is_overlap for randomized_hband = " << difftime(end, start) << endl;
	
	if(pqueue.size() == 0)
		return -1;
		
	max_alignment.clear();
	max_local_match = pqueue.top() / MAXLEN;
	while(!pqueue.empty() && count < MAXTRY)
	{
		dequeue_value = pqueue.top();
		if(DEBUG == 7)
			cout << pqueue.size() << ") dequeue value = " << dequeue_value << endl;
		index = MAXLEN - (dequeue_value % MAXLEN);
		//if(index == MAXLEN) index = 1;
		match = dequeue_value / MAXLEN;
		if(match < max_local_match * HRATIO)
			break;

		if(DEBUG == 6)
			cout << count << ") starting index = " << index << " with matching " << match << endl;
		//update in edit_distance.cpp in github	
		xstring = str1.substr(index - 1, min(str2.length(), str1.length() - index + 1));
		ystring = str2.substr(0, min(xstring.length(), str2.length()));
		//cout << "xstring = " << xstring << endl;
		//cout << "ystring = " << ystring << endl;
		int str1_end, str2_end;
		match = find_kband_similarity(xstring, ystring, alignment, str1_end, str2_end);
		
		if(DEBUG == 5)	printf("The line number is (kband): %d\n", __LINE__);

		if(DEBUG == 6)
		{
			cout << "global alignment score = " << match << endl;	
			//print_alignment(alignment);
		}

		if(global_match <= match && match >= 0)
		{
			global_match = match;
			kband_ratio = 1.0 * match / xstring.length();
			max_alignment.swap(alignment);

			vector<pair<char, char> > vec_seq(max_alignment);
			long long dir_match = match;
			/*
			if(match < 0)
				dir_match = -match * MINLEN + 0;
			else
				dir_match = match * MINLEN + 1;
			*/
			
			dir_match = dir_match * MINLEN + dir;
			dir_match = dir_match * MAXLEN + index;

			assert(dir_match >= 0);
			sequences.push(make_pair(dir_match, vec_seq));
		}
		
		count += 1;	
		alignment.clear();
		pqueue.pop();
	}

	if(kband_ratio < KBAND_THRESHOLD)
		return -1;

	if(FOREACHDIR == 1)
	{
		//cout << "****************************************************************" << endl;
		cout << "maximal alignment with global alignment score = " << global_match << endl;
		print_alignments(max_alignment);
	}

	max_alignment.clear();
	return global_match;
}

