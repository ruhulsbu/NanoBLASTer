#ifndef CONSTANT_H
#define CONSTANT_H

#define SAM_FORMAT 1
#define DEBUG -99		//Debugging code
#define BASE 4			//Number of bases
#define LOCAL 1			//1 for Global alignment; 0 for Local alignment
#define CIRCULAR 0		//0 for circular genome; 1 for non-circular genome


#define FF 1			//Forward forward alignment
#define FR 2			//Forward reverse alignment
#define RF 3			//Reverse forward alignment
#define RR 4			//Reverse reverse alignment


#define GAP_OPEN 2		//Additional penalty for gap opening
#define GAP 6			//Penalty for deletions or insertions (positive)
#define WEIGHT 5		//Reward for matching base
#define MISMATCH -4		//Penalty for mismatch (negative)


#define UP 9			//Symbol for up arrow in alignment path
#define DIAG 6			//Symbol for diagonal arrow in alignment path
#define BACK 3			//symbol for back arrow in alignment path
#define BREAKAT 60		//Number for bases printed to show alignment


#define MATCH 1			//Symbol to calculate matching bases for SAM file
#define INSERT 2		//Symbol to calculate insertions for SAM file 
#define DELETE 3		//Symbol to calculate deletions for SAM file
#define SUBSTITUTE 4		//Symbol to calculate mismatches for SAM file
#define NCHAR 5			//Symbol to calculate N bases in alignment
#define IGNORE 6		//Symbol to ignore the bases
#define MAT_NEXT 2		//Reward for additional matching base 
#define IGNOREN 1		//Ignore the base N


#define MIN_READ_LEN 100	//Minimum read length for alignment
#define MAX_READ_LEN 200000	//Maximum read length for alignment
#define SEED_TUPLE_DIST 10000	//Seed tuple distance in clustering
#define MULTIPLIER 10000000LL	//Used to store x, y as (x * MAXLEN + y)
#define INTERVAL 1


#define PRIM_ANCHOR_WORD 5	//Size of window for primary anchor
#define SECND_ANCHOR_WORD 3	//Size of window for secondary anchor


#define PRIM_ANCHOR_PIDENT 60.0		//Percentage of identity match for primary ANCHOR
#define SECND_ANCHOR_PIDENT 0.65	//Percentage of identity match for secondary ANCHOR


#define FRAGMENT_SIZE 400	//Size of block for sequence alignment
#define KBAND_PERCENT_MATCH 55	//Percentage of identity match for alignment

#define EXPLICIT_SLOPE_90_CHECK 30
#define EXPLICIT_SLOPE_80_CHECK 25
#define EXPLICIT_SLOPE_70_CHECK 20
#define EXPLICIT_SLOPE_50_CHECK 10
#endif
