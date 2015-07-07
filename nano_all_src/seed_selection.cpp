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
                if(DEBUG == 7) 
			cout << "Putting Value = " << lis_vector[i] << endl;
                if(lis_vector[indexing.back()] < lis_vector[i])
                {
                        parent[i] = indexing.back();
                        indexing.push_back(i);
                        if(DEBUG == 7) 
				cout << "Indexing Back = " << indexing.back() << endl;
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

	if(DEBUG == 7)
	{
	        cout << "Indexing Size = " << indexing.size() << endl;
        	for(i = 0; i < indexing.size(); i++)
                	cout << indexing[i] << endl;
	}

	if(DEBUG == 7)
        	cout << "Parent array = " << parent.size() << endl;
        for(u = indexing.size(), v = indexing.back(); u-- ; v = parent[v])
        {
                indexing[u] = v;
		if(DEBUG == 7)
                	cout << u <<") Parent of v = " << v << " is parent[v] = " << parent[v] << endl;
        }

}

int count_kmer_match(string& str1, string& str2, int start, int end)
{
	clock_t t_start, t_end;

	string kmer;
	int hband, i, k, sum;
	vector<int> lis_vector;
	vector<int> indexing;
	//cout << "Length of Str1 = " << str1.length() << ", and Str2 = " << str2.length() << endl;

	t_start = clock();

	hband = min(HBAND, str2.length());

	if(DEBUG == 7)	
	{
		printf("The line number is (inside hband): %d\n", __LINE__);
	}	

/*--------------------------------------------------------------------------------------------------------------------*/

	unordered_map<long, int> map;
	vector<pair<long, int> > vecmap;
	for(int k = 0; k < SEED; k++)
		vecmap.push_back(make_pair(-1, -1));

	bool flag = true;
	long hash_key = 0;
	int map_val; 
	long base_power_kmer = pow(BASE, WINDOW - 1);

	for(int k = 0; k < str2.length(); k++)
	{
		if(flag == true && k + WINDOW > str2.length())
			break;
		for(int l = k, e = k; l < e + WINDOW - 1 && flag == true; l++, k++)
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
			map[hash_key] = k - WINDOW + 1;
		}
		else //now 05-25-15
			map[hash_key] = -1; //now 05-25-15
		//cout << "For character " << read.at(k) << ", at position = " << k <<
		//		", the hash_key = " << hash_key << endl;
	
		vecmap[k] = make_pair(hash_key, k - WINDOW + 1);
		map_val = map_value(str2.at(k - WINDOW + 1));
		hash_key = hash_key - base_power_kmer * map_val;
	}

	flag = true;
	hash_key = 0;
	base_power_kmer = pow(BASE, WINDOW - 1);
	//cout << "LIS Pref for Ref From " << (start - 1) << ", And End = " << (end - 1) << endl;

	for(int k = start - 1, i = 0; k < end - 1; k++, i++)
	{
		if(flag == true && k + WINDOW > end)
			break;
		for(int l = k, e = k; l < e + WINDOW - 1 && flag == true; l++, k++, i++)
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
			cout << "Read Index = " << (map[hash_key]) << ", Ref Index = " << (k - WINDOW + 1) << endl;
			cout << "Read = " << str2.substr(map[hash_key], WINDOW) << endl;
			cout << "Refe = " << str1.substr(k - WINDOW + 1, WINDOW) << endl;
			*/
                        if(map[hash_key] != -1 && (abs(map[hash_key] - (i - WINDOW + 1)) < 5))//5 not 4
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
	
		map_val = map_value(str1.at(k - WINDOW + 1));
		hash_key = hash_key - base_power_kmer * map_val;
	}

/*--------------------------------------------------------------------------------------------------------------------*/

	int lis_size = lis_vector.size();
	t_end = clock();
	t_kmer_count += t_end - t_start;

	t_start = clock();
	if(lis_size < 10) //== 0// < 6 //8 //10
		//return false;
		return -1 * lis_size;

	find_lis_vector(lis_vector, indexing);

	//if(lis_vector.size() < 6)
	//	return false;

	if(DEBUG == 6)
	{
		for (i = 0; i < indexing.size(); i++)
		{
                	printf("\t%d", lis_vector[indexing[i]]);
        		if((i +1) % 12 == 0)
        			printf("\n");
		}
		printf("\n");
	}

	k = lis_vector[indexing[0]];
	sum = WINDOW;
	for(i = 1; i < indexing.size(); i++)
	{
		if(lis_vector[indexing[i]] - k < WINDOW)
		{
			sum = sum + (lis_vector[indexing[i]] - k);
		}
		else
			sum = sum + WINDOW;

		k = lis_vector[indexing[i]];
	}

	cout << "LIS SUM for WINDOW = " << WINDOW << ", and HBAND = " << hband << ", is = " << sum << endl;
	lis_vector.clear();
	indexing.clear();

	t_end = clock();
	t_lis += t_end - t_start;

	if((100.0 * sum / hband)  < PANCHOR_PIDENT)//0.65
	{
		return -1 * sum;//lis_size;
		/*		
		if(lis_size >= 15)
		{
			vector<pair<char, char> > alignment;
        		int str1_end, str2_end, count = 0;
        		string kband_str1 = str1.substr(start + 11, SEED - 11);
        		string kband_str2 = str2.substr(11);
        		find_kband_similarity(kband_str1, kband_str2, alignment, str1_end, str2_end);
        		for(int i =0; i < alignment.size(); i++)
        		{
                		if(alignment[i].first == alignment[i].second)
                		{
                        		count += 1;
                		}
        		}

        		cout << "Maximum number of matching = " << (count + 11) << endl;
        		if((100.0 * (count + 11)) / str2.length() < 71.0)
                		return false;
        		else
                		return true;
		}
		else
		
			return false;
		*/
	}
	else
	{
		return lis_size;
		/*	
		cout << "SUBJECT = " << str1.substr(start - 1, end - start) << endl;
		for(int i = start - 1, k = 0; i < end - WINDOW; i++, k++)
			cout << "\t" << k << " = " << str1.substr(i, WINDOW);
		cout << endl;
		cout << "QUERY   = " << str2 << endl;
		for(int i = 0; i < str2.length() - WINDOW + 1; i++)
			cout << "\t" << i << " = " << str2.substr(i, WINDOW);
		cout << endl;
		
		return true;
		*/
	}
}

int count_kmer_match_prev(string& str1, string& str2, int start, int end)
{
	string kmer;
	int hband, i, k, sum;
	vector<int> lis_vector;
	vector<int> indexing;
	vector<pair<string, int> > vecmap;
	unordered_map<string, int> map;
	//cout << "Length of Str1 = " << str1.length() << ", and Str2 = " << str2.length() << endl;

	hband = min(HBAND, str2.length());

	if(DEBUG == 7)	
	{
		printf("The line number is (inside hband): %d\n", __LINE__);
	}	

	for(i = 0; i < hband - WINDOW + 1; i++)
        {
                kmer = str2.substr(i, WINDOW);
		cout << "kmer  read <" << i << "> = " << kmer << endl;
                if(map.find(kmer) == map.end())
                {
                        map[kmer] = (i);
                }
		else 
			map[kmer] = -1;
		vecmap.push_back(make_pair(kmer, i));
        }

	if(DEBUG == 7)	
	{
		printf("The line number is (inside hband): %d\n", __LINE__);
	}	

	for(i = start - 1, k = 0; i < end - WINDOW; i++, k++)
	{
		kmer = str1.substr(i, WINDOW);	
		cout << "kmer ref  <" << i << "> = " << kmer << endl;
		if(map.find(kmer) != map.end())
		{
			if(map[kmer] != -1 && (abs(map[kmer] - k) < 5))
			{
				lis_vector.push_back(map[kmer]);
			}
			else
			{
				for(int j = max(0, k - 3); j < min(k + 3, vecmap.size()); j++)
				{
					if(vecmap[j].second != -1 && kmer.compare(vecmap[j].first) == 0)
						//kmer.find(vecmap[j].first) != std::string::npos)
					{
						lis_vector.push_back(vecmap[j].second);
						vecmap[j].second = -1;
						break;
					}
				}
			}
		}	
	
	}

	int lis_size = lis_vector.size();
	if(lis_size < 10) // < 8
		return -1 * lis_size;

	find_lis_vector(lis_vector, indexing);

	if(DEBUG == 6)
	{
		for (i = 0; i < indexing.size(); i++)
		{
                	printf("\t%d", lis_vector[indexing[i]]);
        		if((i +1) % 11 == 0)
        			printf("\n");
		}
		printf("\n");
	}

	k = lis_vector[indexing[0]];
	sum = WINDOW;
	for(i = 1; i < indexing.size(); i++)
	{
		if(lis_vector[indexing[i]] - k < WINDOW)
		{
			sum = sum + (lis_vector[indexing[i]] - k);
		}
		else
			sum = sum + WINDOW;

		k = lis_vector[indexing[i]];
	}

	cout << "LIS SUM for WINDOW = " << WINDOW << ", and HBAND = " << hband << ", is = " << sum << endl;
	lis_vector.clear();
	indexing.clear();

	if((100.0 * sum / hband)  < 60.0)
		return -1 * lis_size;
	else
		return lis_size;
}



int find_hband_similarity(string& str1, string& str2, priority_queue<long>& pqueue, int start, int end)
{
	//cout << "Reference Index = " << (start - 1) << endl;
	cout << "similarity of str1 = " << str1.substr(start - 1, SEED) << endl;
	cout << "similarity of str2 = " << str2 << endl << endl;

	/*
	long placement;
	int i, j, k, match, norm_match;
	int hband, index = -1, max_local_match = -1;
	//int matrix[HBAND + 1][str1.length() + 1];
	float hband_ratio = 0.0;
	*/
	if(DEBUG == 7)	
	{
		printf("The line number is (inside hband): %d\n", __LINE__);
	}	

	if(str2.length() * 1 <  SEED)//2
		return -1;
						
	//if(count_kmer_match(str1, str2, start, end) == true)
		
	int return_val = count_kmer_match(str1, str2, start, end);
		
  	if(return_val > 0)
		return 1;
	else
	{
		if(HIGHSITIVE == false)
			return -1;
	}
	//cout << "The LIS Return Value = " << return_val << endl << endl;
	//if(return_val * -1 < 10)
	if(-1.0 * return_val / SEED < 0.50)
		return -1;
		
	//03-26-15
	vector<pair<char, char> > alignment;
	int str1_end, str2_end, count = 0;
	string kband_str1 = str1.substr(start - 1, SEED);
	string kband_str2 = str2;

	//cout << "kband_str1 == " << kband_str1 << endl;
	//cout << "kband_str2 == " << kband_str2 << endl;
	
	///////////////////////////////////////////////////////////////////////////////
	//find_kband_similarity(kband_str1, kband_str2, alignment, str1_end, str2_end);
	int row, column;
        int kband, match, preoff = -1, shift;
        int i, j, h, k, offset;
        float len_ratio;
        string swap;
	/*	
        if(str1.length() < str2.length())
        {
                swap = str1;
                str1 = str2;
                str2 = swap;
        }
	*/

	//cout << "str1 = " << str1 << endl;
	//cout << "str2 = " << str2 << endl;
        kband = min(3, (kband_str2.length()) / 2);
	//kband = min(FRAGMENT_SIZE / 4, str2.length() / 2);//03-26-15//70

        if(DEBUG == 7)
        {
                printf("kband in use = %d\n", kband);
                cout << "First String = " << (kband_str1.length() + 1) << " and Second String = "
                        << (kband_str2.length() + 1) << " and KBAND = " << (2 * kband + 2) << endl;
        }
        row = kband_str1.length() + 1;
        column = 2 * kband + 2;

        if(DEBUG == 7)  printf("The line number is (inside kband): %d\n", __LINE__);
	/*
	cell **matrix = NULL;
        matrix = new cell *[str1.length() + 1];
        for(i = 0; i < str1.length() + 1; i++)
        {
                matrix[i] = new cell[2 * kband + 2];
                memset(matrix[i], 0, (sizeof(cell) * (2 * kband + 2)));
        }
	*/
        if(DEBUG == 7)  printf("The line number is (kband memory dec): %d\n", __LINE__);
	/*
        for(i = 1; i < row; i++)
                for(j = 1; j < column; j++)
                        matrix[i][j].cost = -(i + j) * GAP * LOCAL;
	*/
        len_ratio = 1.0 * kband_str2.length() / kband_str1.length();
        if(DEBUG == 1)
                printf("The Length Ratio = %f\n", len_ratio);
	
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

        if(DEBUG == 7)
        {
                printf("str1.len = %lu and str2.len = %lu\n", kband_str1.length(), kband_str2.length());
                printf("The value of KBAND is = %d\n", (2 * KBAND + 2));
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

                if(DEBUG == 1)
                        cout << "i = " << i << endl;
                if(DEBUG == 1)
                        printf("The offset value is %d\n", offset);

                if(DEBUG == 1)
                        for(k = 1; k < offset - kband; k++)
                                cout << "\t\t";

		start_col = max(1 - offset, -kband);
		end_col = min(kband_str2.length() - offset, kband);

                //for(h = -kband; h <= kband; h++)
                affinity = 0;
                for(h = start_col; h <= end_col; h++)
                {
                        k = offset + h;

                        if(k >= 1 && k <= kband_str2.length())
                        {
                                //match = similarity(str1.at(i - 1), str2.at(k - 1));
				
				if(kband_str1.at(i - 1) == kband_str2.at(k - 1))
				{
					match = WEIGHT;
					
					if(i > 1 && k > 1)//03-26-15
					{
						if(kband_str1.at(i - 2) == kband_str2.at(k - 2))
						{
							match = WEIGHT + MAT_OPEN;
							affinity += 1;
							if(affinity == 6 || affinity == 11)
								match += MAT_OPEN;
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
							match = MISMATCH - MAT_OPEN;
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

                                if(DEBUG == 1)
                                {
					//cout << "@(" << i << "," << j << "," << k << "): " << matrix[i][j] << "\t";
                                        //cout << "@(" << str1.at(i - 1) << ", " << str2.at(k - 1) << "): " << matrix[i][j] << "\t";
                                        cout << "(" << kband_str1.at(i - 1) << "," << kband_str2.at(k - 1) << "):"
                                                << "(" << i << "," << j << "):" << matrix[i][j].cost << "\t";
                                }

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

		//if(max_match + row - i < 0.51 * (max_length + row - i))
		//	break;

                preoff = offset; //to prevent updating same row values
		if(j + 1 < column)
                        matrix[i][j + 1].cost = -(i + j + 1) * GAP * LOCAL;


                if(DEBUG == 1) {
                        cout << "\n";
                        cout << "-------------------------------------------------------------------------------" << endl;
                }


        }

        if(DEBUG == 5)
                cout << "$(" << kband_str1.length() << "," << (kband + 1) << ") Score=" << matrix[kband_str1.length()][kband + 1].cost << endl;
        
	//cout << "Counted max_row = " << max_row << ", and max_column = " << max_col << endl;
	max_row = SEED;
	max_col = kband + 1; 
	//cout << "Final Max_Row = " << max_row << ", and Max_Column = " << max_col << endl;

        //no use now//print_path(path, row - 1, column / 2, str1, str2, alignment);//delete
        {
        	//print_path_back(path, row - 1, column / 2, str1, str2, alignment);
        	//cout << "Total score = " << matrix[max_row][max_col] << ", with kband_score = " << max_score << endl;
        	//cout << "Max score found at the row = " << str1_end << ", and at the col = " << str2_end << endl;
                print_path_cell_back(matrix, max_row, max_col, kband_str1, kband_str2, alignment);
        }


        if(DEBUG == 7)
                printf("The line number is (calling print_path): %d\n", __LINE__);

        if(DEBUG == 1)
                cout << "Printing Path Completed" << endl;

        match = matrix[kband_str1.length()][kband + 1].cost;
	/*
        for(i = 0; i < row; i++)
        {
                delete [] matrix[i];
                matrix[i] = NULL;
        }

        delete [] matrix;

        matrix = NULL;
	*/
        //if(DEBUG == 5) 
        printf("The KBAND matching score = %d\n", match);
        if(DEBUG == 5)  printf("The line number is (end of kband: %d\n", __LINE__);
	/*
        if(OPTIMIZE == 1)
                return optimize_path(alignment);
        else
                return match;
	*/

	///////////////////////////////////////////////////////////////////////////////
	int preoccur = 0;
	int totalScore = 0;
	int windowCount = 0;
	
	count = 0;
	for(int i = 0; i < alignment.size(); i++)
	{
		if(alignment[i].first == alignment[i].second)
		{
			//count += 1;
			preoccur += 1;
			
			if(preoccur == ANCHOR_WORD)
			{
				count += ANCHOR_WORD;
				totalScore += (preoccur * (WEIGHT + MAT_OPEN) - MAT_OPEN);
				//windowCount += 1;
			}
			if(preoccur > ANCHOR_WORD)
			{
				count += 1;
				totalScore += (WEIGHT + MAT_OPEN);
			}
			/*//good for 45 5 55
			if(preoccur == WINDOW)
				totalScore += (WINDOW * (WEIGHT + MAT_OPEN) - MAT_OPEN);
			if(preoccur > WINDOW)
				totalScore += (WEIGHT + MAT_OPEN);
			*/
		}
		else
		{
			preoccur = 0;
			windowCount += 1;
		}
	}
	
	for(int i = 0; i < alignment.size(); i++)
		cout << alignment[i].first;
	cout << endl;
	for(int i = 0; i < alignment.size(); i++)
		cout << alignment[i].second;
        cout << endl;
	cout << "Total kband_score Found = " << match << endl;
	//cout << "Total Window Count = " << windowCount << endl;
	cout << "Total Score accumulated = " << totalScore << endl;
	cout << "Maximum number of matching = " << count << endl << endl;	
	
	
	double quality = 4.5;//4.60;//4.35;//4.50;//4.40;4.6	
	double pident = 0.70;//0.70;//0.65;//0.70;//0.65;70
	
	if(fabs(SANCHOR_PIDENT - 0.70) < 0.3 || SANCHOR_PIDENT > 0.70)
	{
		pident = 0.70;
		quality = 4.5;
	}
	else if(fabs(SANCHOR_PIDENT - 0.65) < 0.3)
	{
		pident = 0.65;
		quality = 4.4;
	}
	else
	{
		pident = SANCHOR_PIDENT;
		quality = 4.35;
	}

	if((1.00 * count  / SEED) < pident || (1.0 * totalScore / SEED) < quality)
		return -1;
	else
		return 1;

}

