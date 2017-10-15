// File general.h - Macros, prototypes and constant of the HMM package
#ifndef _GENERAL_H_
#define _GENERAL_H_


// New types defined 
// typedef unsigned long int word;   // machine word (32 bit on windows - pc). Already defined in 'general.h'

typedef unsigned char byte;       // a byte

//#define ANSI

// Macroes 
#define BIT(block, num)  ( ((block) >> (num)) & (1UL) )  // get the i-th bit from a word


// inverse the order of the bits in a word
#define INVERSE(block) \
	(block) = (((block)&0x0000FFFF) << 16 ) ^ (((block) >> 16)&0x0000FFFF); \
	(block) = (((block)&0x00FF00FF) << 8 ) ^ (((block) >> 8)&0x00FF00FF); \
	(block) = (((block)&0x0F0F0F0F) << 4 ) ^ (((block) >> 4)&0x0F0F0F0F); \
	(block) = (((block)&0x33333333) << 2 ) ^ (((block) >> 2)&0x33333333); \
	(block) = (((block)&0x55555555) << 1 ) ^ (((block) >> 1)&0x55555555); 

// inverse the even and the odd bits in a block
#define INV_PAR(block) \
	(block) = (((block)&0x55555555) << 1 ) ^ (((block) >> 1)&0x55555555); 

// inverse the order of bytes in a word (this is for changing indians)
#define INDY(block) \
	(block) = (((block)&0x0000FFFF) << 16 ) ^ (((block) >> 16)&0x0000FFFF); \
	(block) = (((block)&0x00FF00FF) << 8 ) ^ (((block) >> 8)&0x00FF00FF); 




// swap elements 
#define SWAP(A, B, TEMP)	\
			(TEMP) = (A);   \
			(A) = (B);      \
			(B) = (TEMP);




#define SWAP_S(addr) \
	(addr) = ((addr)&0x20) ^ (((addr)&0x1) << 4) ^ (((addr)&0x1E) >> 1);

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))


#define ABS(x)  ((x) > 0 ?  (x) : (-(x)))


#ifndef TRUE
#define TRUE 0
#endif
#ifndef FALSE
#define FALSE 1
#endif


#define MAX_MODEL_SIZE 6

#define WORD_SIZE 32
#define HALF_WORD_SIZE (WORD_SIZE/2)

#define SIZE_IN_WORDS(A) ((A) >> 4)   /* A usually represents the size in basepairs */

// Get the X'th base from the word A. 0 <= X <= 15
#define GET_BASE(A, X) ( ((A) >> (2*(X)))&0x3 )

// Get the substring containind the base-pairs X to Y : 0 <= X <= Y <= 15
#define GET_BASES(A, X, Y) ( ((A) >> (2*(X)))&0x3 )


// A - 00 T - 01 C - 10 G - 11
// The complement base A<->T and C<->G
#define COMP_BASE(A) ((A) ^ 1) 

// end of line in ASCII
#define EOLN 10

// TAB in ASCII
#define TAB 9

#define SEQ_LEN 100


// Maximal length in words - corresponding to 64*16 = 1024 basepairs
#define MAX_SEQ_LEN 64 

#define MAX_NAME_LENGTH 40

#define MAX_SAMPLES 300

#define RAMAS_SAMPLES 256

// 12980
#define HUMAN13K_GENES 1000  

#define MAX_GENES_NUM HUMAN13K_GENES 

#define EXP_GENES_NUM 7129

#define SLAK_ITERS 50

#define MAX_L 30

#define L_GENES 235

#define L_TF  24




// Maximal number of genes 
#define MAX_NUM_SEQS 1000  

// The gene struct. This contains all the neccessary information needed
// on the specific gene.
typedef struct {


char *seq;        // The sequence of nucleotide. Is this the gene or the promoter ???
long seq_len;     // The length of the sequence 
char *name;       // The name of the gene 
long locus;       // The locus link number 
long NT;          // the accession number NT

int   *seq_exp;   

double *exp_vec;  // Expression Vector for this gene

long exp_len;     // ????

char *symbol;   
char *afy_name;
char *accession;


} gene;


#define NUM_BASEPAIRS_PER_HMM_TRANS 250000.0

#define PRINT_IN_ROW 20

// Number of different values for X  (better take multiplies of 8 !!!)
#define MAX_X_VALS 5 

// Number of different values / Mixtures for Y 
#define MAX_Y_VALS 5  


// Maximal power of M to compute 
#define MAX_POWER_OF_TRANS_MATRIX 50


// Discrete and Continunuos flags
#define DISCRETE 0
#define CONTINUOUS 1

// Number of steps we do until we perform scaling !!
#define DO_SCALE 1
  


#define ONE_OVER_2SQRT_PI  0.3989422802
#define ONE_OVER_2SQRT_PI_LOG (-0.9189385337)

#define EPSILON 0.000001


// This is the number of iterations afterwhich we check if our score is improved !!! 
#define FIRST_ITERS 10
// This number indicates the closeness forwhich we continue to iterate
#define IMPROVEMENT_FRACTION 0.01   

// Maximum possible number of starting points on EM
#define MAX_STARTING_POINTS 50


// This structure contains all the needed information on the HMM model.
// Currently Supported : X - discrete (finite states)
//						 Y - discrete / Mixture of Gaussians
typedef struct {

long x_dim; // dimension of the hidden variable x
long y_dim; // dimension of the observed variable y (in discrete case) or number of mixtures (in continunuos case)

long y_type;   // 0 - Discrete 1 - Continunuos Mixture of Gaussian


double M[MAX_X_VALS][MAX_X_VALS];  // The transition matrix. M[i][j] = Prob (X_n = j | X_{n-1} = i)
double PI[MAX_X_VALS];  // The initial distribution for the i-th

double N[MAX_X_VALS][MAX_Y_VALS];  // The ommision matrix  N[i][j] = Prob (X_n = j | X_{n-1} = i)
double MEW[MAX_X_VALS][MAX_Y_VALS];  // The means of the gaussians
double SIGMA[MAX_X_VALS][MAX_Y_VALS];   // the variances of the gaussians


// Cumulative sums for better performance 
double M_cum[MAX_X_VALS][MAX_X_VALS];  // The transition matrix. M[i][j] = Prob (X_n = j | X_{n-1} = i)
double PI_cum[MAX_X_VALS];  // The initial distribution for the i-th
double N_cum[MAX_X_VALS][MAX_Y_VALS];  // The ommision matrix  N[i][j] = Prob (X_n = j | X_{n-1} = i)


// logs for better performance
double M_log[MAX_X_VALS][MAX_X_VALS];  // The transition matrix. M[i][j] = Prob (X_n = j | X_{n-1} = i)
double PI_log[MAX_X_VALS];  // The initial distribution for the i-th
double N_log[MAX_X_VALS][MAX_Y_VALS];  // The ommision matrix  N[i][j] = Prob (X_n = j | X_{n-1} = i)

double SIGMA_log[MAX_X_VALS][MAX_Y_VALS];  // the standard deviations logs 


double SIGMA_inv[MAX_X_VALS][MAX_Y_VALS];  // the standard deviations inverses

long cum_flag;  // flag saying if to compute cums
long log_flag; // flag saying if to compute logs


// Now the updating parameters flags. Used if some of the parameters are known
long update_M_flag;   // transition matrix updating
long update_N_flag;   // emission matrix updating
long update_MEW_flag;   // mean matrix updating 
long update_SIGMA_flag;   // variance matrix updating 


long fold_change_flag;   // New : a flag saying if the distributions represet
						 // fold change. If so, the MEW parameters are set such that 
						 // the lowest level is no more than one, and the highest level 
						 // is no less then one.  

long place_flag;  // If this flag is on, we use different gaussian for every point


double M_POWERS[MAX_POWER_OF_TRANS_MATRIX][MAX_X_VALS][MAX_X_VALS];  // The transition matrix. M[i][j] = Prob (X_n = j | X_{n-1} = i)
double M_POWERS_cum[MAX_POWER_OF_TRANS_MATRIX][MAX_X_VALS][MAX_X_VALS];  // The transition matrix. M[i][j] = Prob (X_n = j | X_{n-1} = i) cumsum


double *place_gaussian_mew[MAX_X_VALS][MAX_Y_VALS];  // the means of the gaussians
double *place_gaussian_sigma[MAX_X_VALS][MAX_Y_VALS];  // the standard deviations of the gaussians

long miss_data; // Note ! Exists in both structures !!! 


double M_upperbounds[MAX_X_VALS][MAX_X_VALS];  // Bounds on the transition matrix. M[i][j] = Prob (X_n = j | X_{n-1} = i)

long use_bounds;  // Flag saying if to use the bounds on the transition matrix

}  hmm_model;





// This structure contains all the data and related things
typedef struct { 

long miss_data;    // a flag saying if there is a missing data in the Y's - also in hmm_model !!!  

double dont_miss_prob;   // prob. of NOT missing a Y letter (currently simulation assumes i.i.d. places)

// data length 
long seq_len; 

long y_type; // also in hmm_model !!! 

// here is the data 
long *x_vec;
double *y_vec;  
long *y_vec_int; 
long *mix_vec; 
long *loc_vec;
long *loc_diff_vec;



} hmm_data;

#endif