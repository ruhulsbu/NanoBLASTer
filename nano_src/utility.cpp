#include "library.h"
#include "constant.h"
#include "global.h"
#include "structure.h"
#include "functions.h"

int map_value(char ch)
{
	if(ch == 'A')
		return 0;
	if(ch == 'C')
		return 1;
	if(ch == 'G')
		return 2;
	if(ch == 'T')
		return 3;

	return -1;
}

void init_matrix()
{
        matrix = new cell *[2 * FRAGMENTSIZE + 5];
        for(int i = 0; i < 2 * FRAGMENTSIZE + 5; i++)
        {
                matrix[i] = new cell[FRAGMENTSIZE + 5];
        }

}

void remove_matrix()
{
	for(int i = 0; i < 2 * FRAGMENTSIZE + 5; i++)
        {
                delete [] matrix[i];
                matrix[i] = NULL;
        }

        delete [] matrix;

        matrix = NULL;
}

int similarity(char x, char y)
{
	if(x == y)
		return 1 * WEIGHT;
	else 
		return MISMATCH;
}

int min(int x, int y)
{
	if(x < y)
		return x;
	else 
		return y;
}

int max(int x, int y) 
{
	if(x < y)
		return y;
	else
		return x;
}

void upper_case(string& str)
{
	char *array = new char [str.length() + 1];
	strcpy(array, str.c_str());
	for(int i = 0; i < str.length(); i++)
		array[i] = toupper(array[i]);
	
	str.clear();
	str = string(array);
	delete [] array;
}

void reverse_str(string &str)
{
	int i, k;
	char ch;
	
	//cout << "Forward = " << str.substr(0, 80) << endl;
	for(i = 0, k = str.length() - 1; i < k; i++, k--)
	{
		ch = str[i];
		str[i] = str[k];
		str[k] = ch;
	}
	//cout << "Reverse = " << str.substr(str.length() - 80, 80) << endl;

}

string reverse_complement(string& str)
{
	int i, j;
	char x, y;
	char *input = new char [str.length() + 1];
	char *output = new char [str.length() + 1];
	strcpy(input, str.c_str());
	
	for(i = str.length() - 1, j = 0; i >= 0; i--, j++)
	{
		x = input[i];
		if(x == 'A')
			y = 'T';
		else if(x == 'T')
			y = 'A';
		else if(x == 'C')
			y = 'G';
		else if(x == 'G')
			y = 'C';
		else 
			y = 'N';//might also be = x here

		output[j] = y;
		//if(DEBUG == 1) cout << x << " " << y << endl; 	
	}

	output[j] = '\0';
	string rc(output);
	//cout << "Reverse complement(" << str << ") = \n\t\t" << rc << endl;
	
	delete [] input;
	delete [] output;

	return rc;
}

int parseLine(char* line)
{
        int i = strlen(line);
        while (*line < '0' || *line > '9') line++;
        line[i-3] = '\0';
        i = atoi(line);
        return i;
}

int getValue()
{ //Note: this value is in KB!
        FILE* file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];
    

        while (fgets(line, 128, file) != NULL){
            if (strncmp(line, "VmRSS:", 6) == 0){
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
}
