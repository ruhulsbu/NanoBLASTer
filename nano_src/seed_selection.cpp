#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

void find_lis_vector(vector<int>& lis_vector, vector<int>& indexing)
{
        vector<int> parent(lis_vector.size());
        int i, begin, middle, end, u, v;

        if(lis_vector.empty())
                return;
        indexing.push_back(0);

        for(i = 1; i < lis_vector.size(); i++)
        {
                //cout << "Putting Value = " << lis_vector[i] << endl;
                if(lis_vector[indexing.back()] < lis_vector[i])
                {
                        parent[i] = indexing.back();
                        indexing.push_back(i);
                        
			//cout << "Indexing Back = " << indexing.back() << endl;
                        continue;
                }

                begin = 0;
                end = indexing.size() - 1;

                while(begin < end)
                {
                        middle = (begin + end) / 2;
                        if(lis_vector[i] > lis_vector[indexing[middle]])
                                begin = middle + 1;
                        else
                                end = middle;
                }

                if(lis_vector[i] < lis_vector[indexing[begin]])
                {
                        if(begin > 0)
                        {
                                parent[i] = indexing[begin - 1];
                        }
                        indexing[begin] = i;
                }
        }

        for(u = indexing.size(), v = indexing.back(); u-- ; v = parent[v])
        {
                indexing[u] = v;
		//cout << u <<") Parent of v = " << v << " is parent[v] = " << parent[v] << endl;
        }

}

int count_kmer_match(string& str1, string& str2, int start, int end)
{
	string kmer; 
	int i, k, sum;
	vector<int> lis_vector;
	vector<int> indexing;

	unordered_map<long, int> map;
	vector<pair<long, int> > vecmap;

	//cout << "Length of Str1 = " << str1.length() << ", and Str2 = " << str2.length() << endl;
	for(int k = 0; k < ANCHOR; k++)
		vecmap.push_back(make_pair(-1, -1));

	bool flag = true;
	long hash_key = 0;
	int map_val; 
	long base_power_kmer = pow(BASE, PRIMANCHORWORD - 1);

	for(int k = 0; k < str2.length(); k++)
	{
		if(flag == true && k + PRIMANCHORWORD > str2.length())
			break;
		for(int l = k, e = k; l < e + PRIMANCHORWORD - 1 && flag == true; l++, k++)
		{
			map_val = map_value(str2.at(l));
			if(map_val == -1)
				break;

			hash_key = (hash_key << 2) + map_val;
			//cout << "For character " << read.at(k) << ", at position = " << k <<
			//	", the hash_key = " << hash_key << endl;
		}

		map_val = map_value(str2.at(k));
		if(map_val == -1)
		{
			//cout << "Encountered invalid character N ########################" << endl;
			flag = true;
			hash_key = 0;
			continue;
		}
		else
			flag = false;

		hash_key = (hash_key << 2) + map_val;

		if(map.find(hash_key) == map.end())
		{
			map[hash_key] = k - PRIMANCHORWORD + 1;
		}
		else 
			map[hash_key] = -1; //now 05-25-15
		//cout << "For character " << read.at(k) << ", at position = " << k <<
		//		", the hash_key = " << hash_key << endl;
	
		vecmap[k] = make_pair(hash_key, k - PRIMANCHORWORD + 1);
		map_val = map_value(str2.at(k - PRIMANCHORWORD + 1));
		hash_key = hash_key - base_power_kmer * map_val;
	}

	flag = true;
	hash_key = 0;
	base_power_kmer = pow(BASE, PRIMANCHORWORD - 1);
	//cout << "LIS Pref for Ref From " << (start - 1) << ", And End = " << (end - 1) << endl;

	for(int k = start - 1, i = 0; k < end - 1; k++, i++)
	{
		if(flag == true && k + PRIMANCHORWORD > end)
			break;
		for(int l = k, e = k; l < e + PRIMANCHORWORD - 1 && flag == true; l++, k++, i++)
		{
			map_val = map_value(str1.at(l));
			if(map_val == -1)
				break;

			hash_key = (hash_key << 2) + map_val;
			//cout << "For character " << read.at(k) << ", at position = " << k <<
			//	", the hash_key = " << hash_key << endl;
		}

		map_val = map_value(str1.at(k));
		if(map_val == -1)
		{
			//cout << "Encountered invalid character N ########################" << endl;
			flag = true;
			hash_key = 0;
			continue;
		}
		else
			flag = false;

		hash_key = (hash_key << 2) + map_val;

		if(map.find(hash_key) != map.end())
                {
			/*
			cout << "Read Index = " << (map[hash_key]) << ", Ref Index = " << (k - PRIMANCHORWORD + 1) << endl;
			cout << "Read = " << str2.substr(map[hash_key], PRIMANCHORWORD) << endl;
			cout << "Refe = " << str1.substr(k - PRIMANCHORWORD + 1, PRIMANCHORWORD) << endl;
			*/
                        if(map[hash_key] != -1 && (abs(map[hash_key] - (i - PRIMANCHORWORD + 1)) < 5))//5 not 4
                        {
                                lis_vector.push_back(map[hash_key]);
				//cout << "Passed" << endl;
                        }
			else
			{
				for(int j = max(0, i - 3); j < min(i + 3, vecmap.size()); j++)
				{
					if(vecmap[j].first == hash_key)
					{
						lis_vector.push_back(vecmap[j].second);
						vecmap[j].first = -1;
						break;
					}
				}

			}
			//cout << endl;
                }

		//cout << "For character " << read.at(k) << ", at position = " << k <<
		//		", the hash_key = " << hash_key << endl;
	
		map_val = map_value(str1.at(k - PRIMANCHORWORD + 1));
		hash_key = hash_key - base_power_kmer * map_val;
	}


	int lis_size = lis_vector.size();
	
	if(lis_size < 10) 
		return -1 * lis_size;

	find_lis_vector(lis_vector, indexing);
	/*
	if(DEBUG == 99)
	{
		for (i = 0; i < indexing.size(); i++)
		{
                	printf("\t%d", lis_vector[indexing[i]]);
        		if((i +1) % 12 == 0)
        			printf("\n");
		}
		printf("\n");
	}
	*/
	k = lis_vector[indexing[0]];
	sum = PRIMANCHORWORD;
	for(i = 1; i < indexing.size(); i++)
	{
		if(lis_vector[indexing[i]] - k < PRIMANCHORWORD)
		{
			sum = sum + (lis_vector[indexing[i]] - k);
		}
		else
			sum = sum + PRIMANCHORWORD;

		k = lis_vector[indexing[i]];
	}

	//cout << "LIS SUM for PRIMANCHORWORD = " << PRIMANCHORWORD << " is = " << sum << endl;
	lis_vector.clear();
	indexing.clear();

	if((100.0 * sum / ANCHOR)  < PRIMANCHORPIDENT)//0.65
	{
		return -1 * sum;//lis_size;
	}
	else
	{
		return lis_size;
	}
}

int find_hband_similarity(string& str1, string& str2, priority_queue<long>& pqueue, int start, int end)
{
	/*
	if(DEBUG == 99)
	{
		cout << "similarity of str1 = " << str1.substr(start - 1, ANCHOR) << endl;
		cout << "similarity of str2 = " << str2 << endl << endl;
	}
	*/
	if(str2.length() * 1 <  ANCHOR)
		return -1;
		
	int return_val = count_kmer_match(str1, str2, start, end);
		
  	if(return_val > 0)
		return 1;
	else
	{
		if(HIGHSITIVE == false)
			return -1;
	}
	//cout << "The LIS Return Value = " << return_val << endl << endl;
	
	if(-1.0 * return_val / ANCHOR < 0.50)
		return -1;
		
	//03-26-15
	vector<pair<char, char> > alignment;
	int str1_end, str2_end, count = 0;
	string kband_str1 = str1.substr(start - 1, ANCHOR);
	string kband_str2 = str2;

	///////////////////////////////////////////////////////////////////////////////
	int row, column;
        int kband, match, preoff = -1, shift;
        int i, j, h, k, offset;
        float len_ratio;
        string swap;

        kband = min(3, (kband_str2.length()) / 2);
        row = kband_str1.length() + 1;
        column = 2 * kband + 2;
        len_ratio = 1.0 * kband_str2.length() / kband_str1.length();
	
        for(i = 0; i < row; i ++)
        {
                matrix[i][0].cost = i * -GAP * LOCAL;// * 0;
                matrix[i][0].dir = UP;
		matrix[i][0].matrix_col = 0;
		matrix[i][0].str2_index = 0;
		matrix[i][0].match = 0;
		matrix[i][0].length = 0;
        }

        for(j = 1; j < column; j++)
        {
                matrix[0][j].cost = j * -GAP * LOCAL;
                matrix[0][j].dir = BACK;
		matrix[0][j].matrix_col = j - 1;
		matrix[0][j].str2_index = j;
		matrix[0][j].match = 0;
		matrix[0][j].length = 0;
        }

        int max_score = 0.0;
        int max_row, max_col;
	int max_length, max_match;
        max_row = max_col = 1;
	int gap_cost = 0;
	str1_end = str2_end = 0;
	max_length = max_match = 0;
	int start_col, end_col;
	int affinity = 0;

        for(i = 1; i < row; i++)
        {
                offset = (int) (len_ratio * i);
		start_col = max(1 - offset, -kband);
		end_col = min(kband_str2.length() - offset, kband);
                affinity = 0;

                for(h = start_col; h <= end_col; h++)
                {
                        k = offset + h;

                        if(k >= 1 && k <= kband_str2.length())
                        {
	
				if(kband_str1.at(i - 1) == kband_str2.at(k - 1))
				{
					match = WEIGHT;
					
					if(i > 1 && k > 1)//03-26-15
					{
						if(kband_str1.at(i - 2) == kband_str2.at(k - 2))
						{
							match = WEIGHT + MAT_NEXT;
							affinity += 1;
							if(affinity == 6 || affinity == 11)
								match += MAT_NEXT;
						}
					}
				}
				else
				{
					affinity = 0;
					match = MISMATCH;
				
					if(i > 1 && k > 1)
					{
						if(kband_str1.at(i - 2) == kband_str2.at(k - 2))
							match = MISMATCH - MAT_NEXT;
					}
				} 
				
                                if(offset > kband + 1)
                                        j = k - (offset - kband) + 1;
                                else
                                        j = k;

                                if(preoff == offset || offset <= kband + 1)
                                        shift = -1;
                                else
                                        shift = offset - preoff -1;

                                assert(i >= 1 && i < row);
                                assert(j >= 1 && j < 2 * kband + 2);


                                matrix[i][j].dir = DIAG;
				matrix[i][j].matrix_col = (j + shift);
				matrix[i][j].str2_index = k;
				matrix[i][j].length = matrix[i - 1][j + shift].length + 1;
				if(match == WEIGHT)
					matrix[i][j].match = matrix[i - 1][j + shift].match + 1;

                                if(LOCAL == 0)
                                        matrix[i][j].cost = max(0, matrix[i - 1][j + shift].cost + match);
                                else
                                        matrix[i][j].cost = matrix[i - 1][j + shift].cost + match;

                                if(j + shift + 1 <= 2 * kband + 1)
                                {
					if(matrix[i - 1][j + shift + 1].dir == DIAG)
                                   		gap_cost = GAP + GAP_OPEN;
                                	else
                                        	gap_cost = GAP;

                                        if(matrix[i][j].cost < matrix[i - 1][j + shift + 1].cost - gap_cost) 
					{ // <= changed to <
				                matrix[i][j].dir = UP;
						matrix[i][j].matrix_col = (j + shift + 1);
						matrix[i][j].str2_index = k;
                                        	matrix[i][j].cost = matrix[i - 1][j + shift + 1].cost - gap_cost;
						matrix[i][j].length = matrix[i - 1][j + shift + 1].length + 1;
						matrix[i][j].match = matrix[i - 1][j + shift + 1].match;
					}
                                }

                                if(j - 1 >= 1)// && str2.at(k - 1) != 'N')
                                {
					if(matrix[i][j - 1].dir == DIAG)
                                        	gap_cost = GAP + GAP_OPEN;
                                	else
                                        	gap_cost = GAP;

                                        if(matrix[i][j].cost < matrix[i][j - 1].cost - gap_cost)
					{ // <= changed to <
                                                matrix[i][j].dir = BACK;
						matrix[i][j].matrix_col = (j - 1);
						matrix[i][j].str2_index = k;
                                        	matrix[i][j].cost = matrix[i][j - 1].cost - gap_cost;
						matrix[i][j].length = matrix[i][j - 1].length + 1;
						matrix[i][j].match = matrix[i][j - 1].match;
					}
                                }
				/*
                                if(DEBUG == 99)
                                {
					//cout << "@(" << i << "," << j << "," << k << "): " << matrix[i][j] << "\t";
                                        //cout << "@(" << str1.at(i - 1) << ", " << str2.at(k - 1) << "): " << matrix[i][j] << "\t";
                                        cout << "(" << kband_str1.at(i - 1) << "," << kband_str2.at(k - 1) << "):"
                                                << "(" << i << "," << j << "):" << matrix[i][j].cost << "\t";
                                }
				*/
                                if(k >= 0.75 * i && k <= 1.25 * i)//03-26-15
					if(max_score <= matrix[i][j].cost ||
						(max_row < i && max_col < j && matrix[i][j].cost > 0.50 *  matrix[i][j].length))
                                	{		
                                        	max_score = matrix[i][j].cost;
                                        	max_row = i;
                                        	max_col = j;
                                        	str1_end = i;
                                        	str2_end = k;
						max_length = matrix[i][j].length;
						max_match = matrix[i][j].match;
						//cout << "Updating str1_end = " << i << ", and str2_end = " << str2_end <<
						//      ", with ratio = " << (100.0 * k / i) << endl;
					}

                        }

                }

                preoff = offset; //to prevent updating same row values
		if(j + 1 < column)
                        matrix[i][j + 1].cost = -(i + j + 1) * GAP * LOCAL;

		/*
                if(DEBUG == 99) 
		{
                        cout << "\n";
                        cout << "-------------------------------------------------------------------------------" << endl;
                }
		*/

        }

 	max_row = ANCHOR;
	max_col = kband + 1; 

	//print_path_back(path, row - 1, column / 2, str1, str2, alignment);
	//cout << "Total score = " << matrix[max_row][max_col] << ", with kband_score = " << max_score << endl;
	//cout << "Max score found at the row = " << str1_end << ", and at the col = " << str2_end << endl;
	trace_path_cell_back(matrix, max_row, max_col, kband_str1, kband_str2, alignment);
        
        match = matrix[kband_str1.length()][kband + 1].cost;
	//printf("The KBAND matching score = %d\n", match);
	///////////////////////////////////////////////////////////////////////////////

	int preoccur = 0;
	int totalScore = 0;
	int windowCount = 0;	
	count = 0;

	for(int i = 0; i < alignment.size(); i++)
	{
		if(alignment[i].first == alignment[i].second)
		{
			preoccur += 1;
			
			if(preoccur == SECNDANCHORWORD)
			{
				count += SECNDANCHORWORD;
				totalScore += (preoccur * (WEIGHT + MAT_NEXT) - MAT_NEXT);
			}
			if(preoccur > SECNDANCHORWORD)
			{
				count += 1;
				totalScore += (WEIGHT + MAT_NEXT);
			}
		}
		else
		{
			preoccur = 0;
			windowCount += 1;
		}
	}
	/*
	if(DEBUG == 99)
	{
		for(int i = 0; i < alignment.size(); i++)
			cout << alignment[i].first;
		cout << endl;
		for(int i = 0; i < alignment.size(); i++)
			cout << alignment[i].second;
		cout << endl;
		cout << "Total kband_score Found = " << match << endl;
		cout << "Total Score accumulated = " << totalScore << endl;
		cout << "Maximum number of matching = " << count << endl << endl;	
	}
	*/
	double quality = 4.5;//4.60;//4.35;//4.50;//4.40;4.6	
	double pident = 0.70;//0.70;//0.65;//0.70;//0.65;70
	
	if(fabs(SECNDANCHORPIDENT - 0.70) < 0.3 || SECNDANCHORPIDENT > 0.70)
	{
		pident = 0.70;
		quality = 4.5;
	}
	else if(fabs(SECNDANCHORPIDENT - 0.65) < 0.3)
	{
		pident = 0.65;
		quality = 4.4;
	}
	else
	{
		pident = SECNDANCHORPIDENT;
		quality = 4.35;
	}

	if((1.00 * count  / ANCHOR) < pident || (1.0 * totalScore / ANCHOR) < quality)
		return -1;
	else
		return 1;

}

