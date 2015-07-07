#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

//04-21-2015
int KMER = 11;		//Size of KMER. Choose all KMER from Chromosome and Target and hash them. 
int SEED = 40;		//For the Set of KMER in due order find the SEED with 55% match
int MINREAD = 0;	//Start from Read No MINREAD
int MAXREAD = 0;	//End before MAXREAD is reached
int BOUNDARY = 10;	//Rank the alignment and take BOUNDARY option based on max SEED found
int MAX_MATCHED = 55;	//Maximum matching score
int MAX_SCORED = 0;
int ALIGNMENT_CNT = 0;
bool HIGHSITIVE = false;

int t_lookup = 0;
int t_seed = 0;
int t_extend = 0;
int t_chain = 0;
int t_lis = 0;
int t_kmer_count = 0;
long base_power_kmer[21][4];
long error_free_seg[1000];
long error_dist[10];

vector<reference_index> refindex;	//vector of reference
ofstream fp_error_free_seg;
ofstream fp_error_dist;
ofstream fp_csv;
ofstream fp_blastn;
ostringstream logstr;

void input_reference(vector<pair<string, string> >& reference, string& ref_file, string& sam_file);


string input_file, output_file, reference_file, param;

void print_help()
{
	cout << "Use the following options to run NanoBLASTer:" << endl;
	cout << "-p: To specify one of the Parameters: C10, C25, or C50" << endl;
	cout << "-r: To specify the name of Reference file" << endl;
	cout << "-i: To specify the name of Reads file" << endl;
	cout << "-o: To specify the name of Output file" << endl;
	cout << "-k: To specify the size of KMER" << endl;
	cout << "-s: To specify the size of SEED" << endl;
	cout << "-b: To specify the min number of Trials" << endl;
	cout << "-n: To specify the Number of reads to be aligned" << endl;
	cout << "-h, or -?: To print this Help information." << endl << endl;
}

void prepare_input(int argc, char *argv[])
{
	char defaultive[] = "C10";
        char sensitive[] = "C25";
        char highsitive[] = "C50";
     
	int c;
	char *cvalue = NULL;

	reference_file = "";
	input_file = "";
	output_file = "";

	while((c = getopt(argc, argv, "?hp:r:i:o:k:s:b:n:")) != -1)
	{
		switch(c)
		{
			case 'p':
				cvalue = optarg;	
				if(strcmp(cvalue, defaultive) == 0)
                		{
                        		KMER = 13;
                        		SEED = 45;
                        		BOUNDARY = 10;
                		}
                		else if(strcmp(cvalue, sensitive) == 0)
                		{
                        		KMER = 11;
                        		SEED = 40;
                        		BOUNDARY = 25;
                		}
                		else if(strcmp(cvalue, highsitive) == 0)
                		{
                        		KMER = 11;
                        		SEED = 40;
                        		BOUNDARY = 50;
                        		HIGHSITIVE = true;
                		}
                		else
                		{
                        		cout << "Wrong Parameter Given!!! Consult Help." << endl;
					exit(0);
                		}
				break;
			case 'r':
				cvalue = optarg;
				reference_file = string(cvalue);
				break;
			case 'i':
				cvalue = optarg;
				input_file = string(cvalue);
				break;
			case 'o':
				cvalue = optarg;
				output_file = string(cvalue);
				break;
			case 'k':
				KMER = (int) atol(optarg);
				break;
			case 's':
				SEED = (int) atol(optarg);
				break;
			case 'b':
				BOUNDARY = (int) atol(optarg);
				break;
			case 'n':
				MAXREAD = (int) atol(optarg);
				break;
			case 'h':
			case '?':
			default:
				print_help();
				exit(0);
						
		}		
	}

	if(input_file.length() == 0)
	{
		cout << "No Input File is Given. Consult Help." << endl;
		exit(0);
	}
	if(reference_file.length() == 0)
	{
		cout << "No Reference File is Given. Consult Help." << endl;
		exit(0);
	}
}

int main(int argc, char *argv[])
{
	/*
	char defaultive[] = "C10";
        char sensitive[] = "C25";
        char highsitive[] = "C50";
     
	if(argc == 2)
        {
                if(strcmp(argv[1], defaultive) == 0)
                {
                        KMER = 13;
                        SEED = 45;
                        BOUNDARY = 10;
                }
                else if(strcmp(argv[1], sensitive) == 0)
                {
                        KMER = 11;
                        SEED = 40;
                        BOUNDARY = 25;
                }
                else if(strcmp(argv[1], highsitive) == 0)
                {
                        KMER = 11;
                        SEED = 40;
                        BOUNDARY = 50;
                        HIGHSITIVE = true;
                }
                else
                {
                        cout << "Wrong Parameter Given!!!" << endl;
                }
        }
	if(argc == 3)
	{
		KMER = (int) atol(argv[1]);
		SEED = (int) atol(argv[2]);
	}
	if(argc == 4)
	{
		KMER = (int) atol(argv[1]);
		SEED = (int) atol(argv[2]);
		BOUNDARY = (int) atol(argv[3]);
	}
	if(argc == 5)
	{
		KMER = (int) atol(argv[1]);
		SEED = (int) atol(argv[2]);
		BOUNDARY = (int) atol(argv[3]);
		MAXREAD = (int) atol(argv[4]);
	}
	*/

	prepare_input(argc, argv);
	cout << "KMER = " << KMER << ", SEED = " << SEED << ", BOUNDARY = " << BOUNDARY << ", MAXREAD = " << MAXREAD << endl;
	cout << "Reference = " << reference_file << ", Input = " << input_file << ", Output = " << output_file << endl;

	clock_t t_start, t_end;
	time_t start, end;
	float t_index;
	int input, cases = 0, overlap, index_time;

	string str1, str2;
	time(&start);
	
	init_matrix();
	if(output_file.length() == 0)
	
		logstr << "test";
	else
		logstr << output_file;

	string log = "detail_alignment-";
	//logstr << SEED << "-" << MINREAD << "-" << MAXREAD;
	log = log + logstr.str() + ".txt";

	log = "/dev/null";
	freopen(log.c_str(), "w", stdout);
	//freopen("ont_detail_alignment.txt", "w", stdout);
	
	//ofstream file("/dev/null");
	//cout.rdbuf(file.rdbuf());

	string csv = "stat-" + logstr.str() + ".csv";
	fp_csv.open(csv.c_str(), ofstream::out);

	string blastn = "output-" + logstr.str() + ".txt";
	fp_blastn.open(blastn.c_str(), ofstream::out);

	string ref_file = reference_file;//"../ref.fa";
	//string ref_file = "testref.fa";
	
	string sam_file = ""; //test-	//4. index
	sam_file = sam_file + logstr.str() + ".sam";

	string read_file = input_file;//"../reads.fa";
	//string read_file = "testreads.fa";
	vector<pair<string, string> > reference;

	string file_error_free_seg = "error-free-seg-" + logstr.str() + ".csv";
	memset(error_free_seg, 0, sizeof(long) * 1000);

	string file_error_dist = "error-dist-" + logstr.str() + ".csv";
	memset(error_dist, 0, sizeof(long) * 1 * 10);
	fp_error_dist.open(file_error_dist.c_str(), ofstream::out);

	input_reference(reference, ref_file, sam_file);
	//reference_index refind[reference.size()];
	
	//long base_power_kmer = pow(BASE, KMER - 1);
	for(int i = 0; i < KMER; i++)
		for(int k = 0; k < BASE; k++ )
			base_power_kmer[i][k] = pow(BASE, i) * k;

	clock_t t_start_str;
	clock_t t_end_str;
	int t_total_str;
	long kmer_choices = (long) (pow(4, KMER) + 9);	
	for(int i = 0; i < reference.size(); i++)
	{
		t_start_str = clock();
		//cout << reference.at(i).first << ": " << reference.at(i).second.size() << endl;
		reference_index refinfo;
		refinfo.name = reference.at(i).first;
		refinfo.ref = reference.at(i).second;
		refinfo.rev = reverse_complement(reference.at(i).second);

		refinfo.index = new int [kmer_choices];//[67108869];//[16777219];// [4194305];
		memset(refinfo.index, -1, (sizeof(int) * kmer_choices));//67108869));//16777219));// 4194305));
		refinfo.revind = new int [kmer_choices];
		memset(refinfo.revind, -1, (sizeof(int) * kmer_choices));

		//t_end_str = clock();
		cout << "ref: " << refinfo.name << " and name = " << refinfo.ref.size() << endl;

		//t_start_str = clock();
		long hash_key = 0;
		int map_val = 0;
		for(int k = 0, i = KMER - 2; k < KMER - 1; k++, i--)
		{
			map_val = map_value(refinfo.ref.at(k));
			hash_key += base_power_kmer[i][map_val];
			
			//cout << "Init character " << refinfo.ref.at(k) << ", at position = " << k <<
			//	", the hash_key = " << hash_key << endl;
		}
		t_end_str = clock();
		t_total_str += t_end_str - t_start_str;

		for(int k = KMER - 1; k < refinfo.ref.length(); k++)
		{
			t_start_str = clock();
			map_val = map_value(refinfo.ref.at(k));
			hash_key = (hash_key << 2) + map_val;
			t_end_str = clock();
			t_total_str += t_end_str - t_start_str;	
		
			//time(&start);
			t_start = clock();
			
			if(refinfo.index[hash_key] == -1)
			{
				vector<int> pos;
				refinfo.position.push_back(pos);
				refinfo.index[hash_key] = refinfo.position.size() - 1;
			} 
					
			refinfo.position[refinfo.index[hash_key]].push_back(k - KMER + 1);

			t_end = clock();
			t_index += t_end - t_start;
		
			map_val = map_value(refinfo.ref.at(k - KMER + 1));
			hash_key -= base_power_kmer[KMER - 1][map_val];
		}

		/////////////////////////////////////////////////////////////////////////////////////

		//t_start_str = clock();
		hash_key = 0;
		map_val = 0;
		for(int k = 0, i = KMER - 2; k < KMER - 1; k++, i--)
		{
			map_val = map_value(refinfo.rev.at(k));
			hash_key += base_power_kmer[i][map_val];
			
			//cout << "Init character " << refinfo.ref.at(k) << ", at position = " << k <<
			//	", the hash_key = " << hash_key << endl;
		}
		t_end_str = clock();
		t_total_str += t_end_str - t_start_str;

		for(int k = KMER - 1; k < refinfo.rev.length(); k++)
		{
			t_start_str = clock();
			map_val = map_value(refinfo.rev.at(k));
			hash_key = (hash_key << 2) + map_val;
			t_end_str = clock();
			t_total_str += t_end_str - t_start_str;	
		
			//time(&start);
			t_start = clock();
			
			if(refinfo.revind[hash_key] == -1)
			{
				vector<int> pos;
				refinfo.position.push_back(pos);
				refinfo.revind[hash_key] = refinfo.position.size() - 1;
			} 
					
			refinfo.position[refinfo.revind[hash_key]].push_back(k - KMER + 1);

			t_end = clock();
			t_index += t_end - t_start;
		
			map_val = map_value(refinfo.rev.at(k - KMER + 1));
			hash_key -= base_power_kmer[KMER - 1][map_val];
		}
	
		cout << "Total index size = " << refinfo.position.size() << endl;	
		refindex.push_back(refinfo);
		t_end_str = clock();	
		t_total_str += t_end_str - t_start_str;
		//cout << "ref: " << refind[i].name << " and name = " << refind[i].ref.size() << endl;
		
		//break;
		
	}

	time(&end);
	index_time = difftime(end, start);

	time(&start);
	
	align_reads(reference, read_file, sam_file, refindex);
	
	time(&end);

	cout << "Total Execution time = " << difftime(end, start) << ", and Index time = " << index_time << endl;// << endl;
	cout << "Total Hashing Time = " << (t_index / CLOCKS_PER_SEC) << endl;
	cout << "Total String Time = " << ((float)t_total_str / CLOCKS_PER_SEC) << endl;
	cout << "Total Lookup Time = " << ((float)t_lookup / CLOCKS_PER_SEC) << endl;
	cout << "Total Seeding Time = " << ((float)t_seed / CLOCKS_PER_SEC) << endl;
	cout << "Total Extending Time = " << ((float)t_extend / CLOCKS_PER_SEC) << endl;
	cout << "Total Chaining Time = " << ((float)t_chain / CLOCKS_PER_SEC) << endl;
	cout << "Total LIS Time = " << ((float)t_lis / CLOCKS_PER_SEC) << endl;
	cout << "Total KMER Count Time = " << ((float)t_kmer_count / CLOCKS_PER_SEC) << endl;

	for(int i = 0; i < reference.size(); i++)
	{
		reference[i].first.clear();
		reference[i].second.clear();
	}

	for(int i = 0; i < refindex.size(); i++)
	{
		refindex[i].ref.clear();
		refindex[i].rev.clear();
		refindex[i].name.clear();
		//refindex[i].index.clear();
		delete [] refindex[i].index;
		delete [] refindex[i].revind;
	}

	remove_matrix();
	fp_blastn.close();
	fp_csv.close();

	fp_error_free_seg.open(file_error_free_seg.c_str(), ofstream::out);
	for(int i = 1; i < 250; i++)
	{
		fp_error_free_seg << i << "," << error_free_seg[i] << endl;
	}
	fp_error_free_seg.close();
	
	/*
	fp_error_dist.open(file_error_dist.c_str(), ofstream::out);
	for(int i = 0; i < 100001; i++)
	{
		fp_error_dist << i << "," << error_dist[i][0] << "," << error_dist[i][1] << 
					"," << error_dist[i][2] << "," << error_dist[i][3] << 
					"," << error_dist[i][4] << "," << error_dist[i][5] <<
					"," << error_dist[i][6] << "," << error_dist[i][7] << 
					"," << error_dist[i][8] << "," << error_dist[i][9] << endl;
		//max_error_seg, total_match, total_insert, total_delete
		//total_substitute, total_nchar, total_ignore, alignment_length, 
		//read_length, reference_length 
					
	}
	*/
	fp_error_dist.close();

	exit(0);
	return 0;
}

void input_reference(vector<pair<string, string> >& reference, string& ref_file, string& sam_file)
{
	ifstream fp_ref;
	ofstream fp_sam;
	
	char *ref = new char[ref_file.length() + 1];
	strcpy(ref, ref_file.c_str());

	char *sam = new char[sam_file.length() + 1];
	strcpy(sam, sam_file.c_str());
	

	fp_ref.open(ref, ifstream::in);
	fp_sam.open(sam, ofstream::out);

	string space = " ";
	string input, output;

	string ref_name;
	string sequence;

	getline(fp_ref, input);

	while(!fp_ref.eof())
	{
		//getline(fp_ref, input);
		size_t find = input.find(space);
		if(find != string::npos)
			ref_name = input.substr(1, find);
		else
			ref_name = input.substr(1);
		sequence = "";		

		while(getline(fp_ref, input))
		{
			if(input.at(0) == '>')
			{
				break;
			}
			sequence += input;
		}

		//if(ref_name.find("chr13") == string::npos)
		//	continue;//added on 03-13-15
		
		upper_case(sequence);		
		reference.push_back(make_pair(ref_name, sequence));
		output = "@SQ\tSN:" + ref_name + "\tLN:";
		fp_sam << output << sequence.length() << endl;

		//break; //added on 03-11-15		
	}

	fp_ref.close();
	fp_sam.close();

	delete [] ref;
	delete [] sam;
}
