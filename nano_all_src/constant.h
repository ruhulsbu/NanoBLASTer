//constant variables

#ifndef CONSTANT_H
#define CONSTANT_H

#define WINDOW 5	//size of window//should be 6
#define ANCHOR_WORD 3
#define KBAND 100	//size of band for string comparisons
#define HBAND 200	//currently unused

#define MAXLEN 10000000LL	//used to store x, y as (x * MAXLEN + y) as long long
#define MINLEN 100		//used to store x, y as (x * MINLEN + y) as long
#define MAXTRY 10		//number of trials to select reference for a target
#define MAXREPEATLEN 1000	//not implemented yet
#define SEEDTUPLEDIST 3000	//not implemented yet
#define SPECIFICITY 0.5		//not implemented yet
#define BOUNDED 0		//used for circular genome

#define GAP 6		//delete or insert	//positive
#define MAT_OPEN 2	//start matching bases
#define GAP_OPEN 2	//start of gap
#define WEIGHT 5	//match bases
#define MISMATCH -4	//substitute		//negative
#define DEBUG 6		//debug code
#define IGNOREN 1

#define UP 9		//symbol for up arrow
#define DIAG 6		//symbol for diagonal arrow
#define BACK 3		//symbol for back arrow
#define BREAKAT 60	//number for printing the string alignment

#define MATCH 1
#define INSERT 2
#define DELETE 3
#define SUBSTITUTE 4
#define NCHAR 5
#define IGNORE 6

#define FF 1		//forward forward alignment
#define FR 2		//forward reverse alignment
#define RF 3		//reverse forward alignment
#define RR 4		//reverse reverse alignment

#define HBAND_THRESHOLD 0.0	
#define KBAND_THRESHOLD 0.0

#define LOCAL 1
#define OPTIMIZE 0
#define FOREACHDIR 0 

#define HRATIO 0.95
#define ERROR 0.5
//#define KMER 11
//#define SEED 45 
#define FULLREAD 1
#define BASE 4

#define MAXREADLEN 200000
#define MINREADLEN 100

#define FRAGMENT_SIZE 400
#define GAP_PERCENT_MATCH 55
#define KBAND_PERCENT_MATCH 55
#define ORIGINAL_PERCENT_MATCH 53
#define CHAIN_PERCENT_MATCH 55
#define ALIGNMENT_COVERAGE 4.99// 4.99f


#define PANCHOR_PIDENT 65.0
#define SANCHOR_PIDENT 0.70

#endif
