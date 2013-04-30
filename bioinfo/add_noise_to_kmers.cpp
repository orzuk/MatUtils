#include <stdlib.h>
#include <stdio.h>
//#include <iostream.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <matrix.h>
#include <limits.h>

#include "matrix.h"
#include "mex.h"
#include "general.h"
#include "hmm_chrom_funcs.h"
#include "dna_utils_new.h"
#include "markov.h"
#include "uthash.h" // New! hash table (is this causing us all our memory problems???) 

#undef uthash_expand_fyi
#define uthash_expand_fyi(tbl) printf("expanded to %d buckets\n", tbl->num_buckets)

// Hash table structures 
struct SeqData {
	int seqpos;
	word seq[MAX_L]; /* need to set the right length */
	UT_hash_handle hh; /* makes this structure hashable */
	int rand_int; /* temp for debug */
};
struct SeqData *hTable=NULL;



///////////////////////////////////////////////////////////////////////////////////////
/// Divide a set of sequences to kmers 
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{

	//	printf("REALLY Start mex function\n");	
	//	fflush(stdout);

	int nRows, nCols, i, j, k;

	long L; 
	word **seqs;
	long *seqs_lens;   // New : There may be different sequences lengths

	double *in, *inn, *out;
	word *in_w, *out_w;
	long kmer_length_in_words; 
	long unique_flag=0; // do we want unique kmers (and keep their indices) 
	long hash_flag; //  = 1; 
	long half_word_size = 4*sizeof(word); // half word size in bits (=word size in nucleotides)
	long word_size = 2*half_word_size;
	long num_coordinates; // counts number of output coordinates 
	long *kmer_position_inds;
	long *kmer_seq_inds; // = new long[total_num_kmers]; 
	long input_position_flag=0; // flag saying 	
	long debug_prints=0; 
	double noise_table[MAX_L][4][4]; 

	/**
	printf("Start mex function\n");	
	fflush(stdout);
	printf("Size of Word=%ld (bits), Size of long=%ld\n", 8*sizeof(word), 8*sizeof(long)); // print size in bits
	fflush(stdout);
	/**/

	/*************************************************  Read Input **********************************************************************/
	/* Check for proper number of arguments. Should be five 5 (last one is double flag) */
	/* Input variables are: L, num_kmers, kmers, noise_table */
	if(nrhs != 4) {
		mexErrMsgTxt("Usage: [L], [num_kmers], [kmers], [noise_table]");
	}

	in = mxGetPr(prhs[0]);   // Get the kmer length (L)  	  
	L	= in[0]; 
	kmer_length_in_words = ceil(double(L)/double(half_word_size)); 
	in = mxGetPr(prhs[1]);   // Get the number of kmers     	  
	long num_kmers = in[0]; 
	long num_words_in_kmer = ceil(double(L) / double(half_word_size)); 
//	printf("Reading kmers and noise table"); 
//	fflush(stdout);


	in_w = (word *)(mxGetData(prhs[2]));   // Get the kmers (in packed form) 	  
	word **kmers = new word*[num_kmers]; // read input sequences 
	for(i = 0; i < num_kmers; i++)
	{
		kmers[i] = new word[num_words_in_kmer];   // allocate memory for every sequence
		for(j = 0; j < num_words_in_kmer; j++)
			kmers[i][j] = in_w[j*num_kmers + i];  // Copy kmers (they should be packed for now !!!)
	}	

	if(debug_prints)
	for(i = 0; i < num_kmers; i++) // loop on all sequences 
		{
			printf("Kmer %ld Is: \n", i);
			print_packed_dna(kmers[i], L, 1); 
		}


	in = mxGetPr(prhs[3]);   // Get the noise table (we use 16xL double matrix)     	  
	for(i = 0; i < L; i++)
		for(j = 0; j < 4; j++)
			for(k = 0; k < 4; k++)
				noise_table[i][j][k] = in[   j*4*L+k*L+i];

		if(debug_prints)
		{
			printf("Noise Table\n"); 
			for(i=0; i<L; i++) // compute cumulative sum 
			{
				printf("I=%ld: ", i); 
				for(j=0; j<4; j++)
				{
					for(k=0; k<4; k++)
						printf("%lf, ", noise_table[i][j][k]); 
					printf("  |  "); 
				}
				printf("\n"); 
			}
		}




	//	if ( mxIsUint32(prhs[0]) != 1)
	//      mexErrMsgTxt("Input must be a 32-bit integer!!! ");




	printf("Obtaining coordinates\n"); 
	fflush(stdout);


	if(debug_prints)
	{	
		printf("Allocate seqs lengths  %ld\n", num_kmers);
		fflush(stdout);
	}




	//	printf("Inputs: kmer-length=%ld, unique-flag=%ld, hash-flag=%ld\n", L, unique_flag, hash_flag); 
	//	fflush(stdout); 



	//	printf("Allocating %ld\n", total_num_kmers);
	//	fflush(stdout);
	word **noisy_kmers = new word*[num_kmers]; 	//  	word **kmers = (word **) mxMalloc( total_num_kmers * sizeof(double));
	for(i=0; i<num_kmers; i++)
		noisy_kmers[i] = new word[num_words_in_kmer];   // allocate memory for every sequence


	//	printf("Allocated main array kmers %ld\n", total_num_kmers);
	//	fflush(stdout);

	//	printf("Allocated main array indices %ld\n", total_num_kmers);
	//	fflush(stdout);


	// print to see that input came right
	if(debug_prints)
	{
		printf("Print input DNA: First kmer=%ld out of %ld\n", kmers[0][0], num_kmers); 
		for(i=0; i<num_kmers; i++)
			print_packed_dna(kmers[i], L, 1);
		fflush(stdout);	
	}
	/*************************************************  Finished Read Input ***************************************************************/


	//	printf("kmer length in words =%ld\n", kmer_length_in_words);



	//	*intersect_start_pos = (double *) mxMalloc( (n1+n2) * sizeof(double));
	/*
	kmers = new word*[total_num_kmers]; // this may be too much - what to do? allocate on the fly? 
	for(i = 0; i < total_num_kmers; i++)
	kmers[i] = new word[kmer_length_in_words]; 
	*/
	//	printf("L is %d, num_seqs %d, seqs_len %d, total_kmers %d kmer-in-word %d unique_flag %ld\n", 
	//		L, num_seqs, seqs_len, total_num_kmers, kmer_length_in_words, unique_flag);


	//	printf("Extract sub-kmers\n");
	//	fflush(stdout);

	//	if(input_position_flag)
	//		return;
	//	else

	add_noise_to_kmers(L, num_kmers, kmers, noise_table, 
		noisy_kmers);





	/*
	for(i=0; i<num_seqs; i++) 
	for(j=0; j<5; j++)
	printf("original word = %lu\n", seqs[i][j]);
	*/


	/*
	for(i=0; i<total_num_kmers; i++) 
	for(j=0; j<kmer_length_in_words; j++)
	printf("%ld\n", kmers[i][j]);
	*/

	//	mwSize mw_dims[2];
	//	mw_dims[0] = mw_total_num_kmers; mw_dims[1] = mw_kmer_length_in_words;



	long return_double = 0;
	if(return_double) // return kmers in doubles
	{
		plhs[0] = mxCreateDoubleMatrix(num_kmers, kmer_length_in_words, mxREAL); //	  plhs[0] = mxCreateNumericArray(2, mw_dims, mxUINT32_CLASS, mxREAL);
		out = mxGetPr(plhs[0]); // output the scores array 
		for(i=0; i<num_kmers; i++) 
			for(j=0; j<kmer_length_in_words; j++)
				out[i+j*num_kmers] = double(noisy_kmers[i][j]); // copy kmers array to output // 0*j*total_num_kmers+
	}
	else // here return long integers
	{
		if(word_size == 32)
			plhs[0] = mxCreateNumericMatrix(num_kmers, kmer_length_in_words, mxUINT32_CLASS, mxREAL);
		else
			plhs[0] = mxCreateNumericMatrix(num_kmers, kmer_length_in_words, mxUINT64_CLASS, mxREAL);

		out_w = (word *)(mxGetData(plhs[0]));
		for(i=0; i<num_kmers; i++) 
		{
			//			printf("Copy kmer %ld\n", i); 
			for(j=0; j<kmer_length_in_words; j++)
				out_w[i+j*num_kmers] = noisy_kmers[i][j]; // copy kmers array to output // 0*j*total_num_kmers+
		}
	}


	for(i = 0; i < num_kmers; i++)
	{
		delete kmers[i]; 			  
		delete noisy_kmers[i];
	}
	delete kmers; //   mxFree(kmers);
	delete noisy_kmers;


}
//#endif // MEX_TO_MATLAB


// free hash memory
void delete_all() {
	struct SeqData *current_user, *tmp;

	HASH_ITER(hh, hTable, current_user, tmp) {
		HASH_DEL(hTable,current_user);  /* delete; users advances to next */
		free(current_user);            /* optional- if you want to free  */
	}
}



// 
// 
////////////////////////////////////////////
//
//  function : add_noise_to_kmers
// 
//  input    : L - length of kmers (in nucleotides)
//			   num_kmers - number of kmers 
//			   kmers - array of kmers (in packed form 2-dim array)
//			   noise_table	- table of size Lx4x4 of noise values. noise_table[i][k][k] = Pr(z_i=k | x_i=j)
//
//  output   :  noisy_kmers - array of noisy kmers 
//
//  purpose  : Add noise to kmers - perform substitutions randomly according to an input noise model (insertions and deletions not supported)			   
//
////////////////////////////////////////////

long add_noise_to_kmers(long L, long num_kmers, word *kmers[MAX_L], double noise_table[MAX_L][4][4], 
						word *noisy_kmers[MAX_L])

{
	long i, j, k;
	long half_word_size = 4*sizeof(word);
	long word_size = 2*half_word_size;

	long num_words_in_kmer = ceil(double(L) / double(half_word_size)); 
	long num_bytes_in_kmer = ceil(double(L) / 4.0); // temp - a bit wasteful but easier to make sure it's correct 
	long num_bytes_in_kmer_ceil = 4 * ceil(double(L) / 16.0); // temp - a bit wasteful but easier to make sure it's correct 
	word cur_nuc;
	long last_word_in_kmer = (L-1)%half_word_size+1; 
	word current_kmer[MAX_L]; // save one kmer 
	long debug_prints=0; 
	long add_int = 0; 
	word ones_mask;
	double cumulative_noise_table[MAX_L][4][4]; 

	double r; // random number 
	if(word_size == 32)
		ones_mask = 0xFFFFFFFF;
	else // assume it's 64 
		ones_mask = 0xFFFFFFFFFFFFFFFF;


	if(debug_prints)
	{
		printf("Start adding noise\n"); 
		fflush(stdout);
	}

	for(i=0; i<L; i++) // compute cumulative sum 
		for(j=0; j<4; j++)
		{
			cumulative_noise_table[i][j][0] = noise_table[i][j][0];
			for(k=1; k<4; k++)
				cumulative_noise_table[i][j][k] = cumulative_noise_table[i][j][k-1] + noise_table[i][j][k];
		}

		if(debug_prints)
		{
			printf("Cumulative Noise MAtrix\n"); 
			for(i=0; i<L; i++) // compute cumulative sum 
			{
				printf("I=%ld: ", i); 
				for(j=0; j<4; j++)
				{
					for(k=0; k<4; k++)
						printf("%lf, ", cumulative_noise_table[i][j][k]); 
					printf("  |  "); 
				}
				printf("\n"); 
			}
		}

		for(i = 0; i < num_kmers; i++) // loop on all sequences 
		{
//			printf("Adding noise to kmer %ld\n", i); 
//			fflush(stdout);
			memcpy(noisy_kmers[i], kmers[i], num_bytes_in_kmer_ceil); // first copy without noise
			for(j = 0; j < L; j++) // loop on all positions 
			{
				cur_nuc = BITS(kmers[i][j/half_word_size], 2*(j%half_word_size), 2*(j%half_word_size)+1); // get kmer represented as two bits
//				printf("---Adding noise to position %ld\n. Cur nuc=%ld. ", j, cur_nuc); 
/*				printf("Noise Matrix used: [%lf %lf %lf %lf] Cumulative:  [%lf %lf %lf %lf]\n", 
					noise_table[j][cur_nuc][0], noise_table[j][cur_nuc][1], 
					noise_table[j][cur_nuc][2], noise_table[j][cur_nuc][3],
					cumulative_noise_table[j][cur_nuc][0], cumulative_noise_table[j][cur_nuc][1], 
					cumulative_noise_table[j][cur_nuc][2], cumulative_noise_table[j][cur_nuc][3]);
				fflush(stdout);
*/
				r = double(rand()) / RAND_MAX; // randomize in [0,1]
//				printf(" Randomized r=%lf\n", r); 
				for(k=0; k<4; k++)
				{
//					printf("------Try nucleotide %ld ", k); 
//					fflush(stdout);
					if(r <= cumulative_noise_table[j][cur_nuc][k]) // passed 
					{
//						printf("nuc[%ld][%ld] %ld -> %ld\n", i, j, cur_nuc, k); 
						if(k != cur_nuc) // noisy kmer different from true kmer
						{
							noisy_kmers[i][j/half_word_size] &= (ones_mask ^ (3UL << (2*(j%half_word_size)))); // get rid of previous nucleotide
							noisy_kmers[i][j/half_word_size] ^= (k << (2*(j%half_word_size))); 
						}
						break; 
					}


				}
			}
		}

		if(debug_prints)
		for(i = 0; i < num_kmers; i++) // loop on all sequences 
		{
			printf("Kmer %ld: [pure,noisy]\n", i);
			print_packed_dna(kmers[i], L, 1); 
			print_packed_dna(noisy_kmers[i], L, 1); 
		}

		return 0; 
} // end of function extract_sub_kmers

////////////// Temp /////////////////////
// len is in nucleotides 
// print_mode: 1 - letters, 0 - just numbers 
void print_packed_dna(word *seqs, long len, long print_mode)
{
	long i, j;
	word tmp;
	long half_word_size = 4*sizeof(word);
	word ones_mask;
	if(half_word_size == 32)
		ones_mask = 0x1F; 
	else
		ones_mask=0xF;


	for(i = 0; i < len; i++)
	{
		//		printf("\nSeq[i/16]=%ld\n", seqs[i/16]);
		tmp = seqs[i/half_word_size];

		tmp = tmp >> (2*(i&ones_mask));
		tmp = tmp&0x3;
		tmp++;

		//		tmp = ((seqs[i/16] << (2* (i&0xF)))&0x3) + 1;
		//		tmp = ((seqs[i/16] >> ((i&0xF)*2))&0x3) + 1;
		if(print_mode == 0)
			printf("%lx ", tmp);
		else
		{
			//			printf("%lx ", tmp);
			switch (tmp)
			{
			case 1:
				printf("A"); break;
			case 2:
				printf("C"); break;
			case 3: 
				printf("G"); break;
			case 4:
				printf("T"); break;
			}
		}
	}
	printf("\n");

}
