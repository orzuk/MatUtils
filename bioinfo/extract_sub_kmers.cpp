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

	int nRows, nCols, i, j, K;

	long L; 
	word **seqs;
	long *seqs_lens;   // New : There may be different sequences lengths
	long num_seqs; 
	//  long seqs_len; 
	//	word **kmers; // output kmers 

	double *in, *inn, *out;
	word *in_w, *out_w;
	long total_num_kmers = 0; 
	long kmer_length_in_words; 
	long unique_flag=0; // do we want unique kmers (and keep their indices) 
	long num_unique_kmers, num_sparse_matrix_elements;
	long hash_flag; //  = 1; 
	long half_word_size = 4*sizeof(word); // half word size in bits (=word size in nucleotides)
	long word_size = 2*half_word_size;
	long num_coordinates; // counts number of output coordinates 
	long *kmer_position_inds;
	long *kmer_seq_inds; // = new long[total_num_kmers]; 
	long input_position_flag=0; // flag saying 	
	long debug_prints=0; 

	/**
	printf("Start mex function\n");	
	fflush(stdout);
	printf("Size of Word=%ld (bits), Size of long=%ld\n", 8*sizeof(word), 8*sizeof(long)); // print size in bits
	fflush(stdout);
	**/

	/*************************************************  Read Input **********************************************************************/
	/* Check for proper number of arguments. Should be five 5 (last one is double flag) */
	/* Input variables are: pwms, seqs, seqs_lens, -10, is_double */
	if(nrhs < 5) {
		mexErrMsgTxt("Usage: [seqs_packed], [lens], [K], [unique-flag], [hash-flag]");
	}


	//		(word *)(mxGetPr(prhs[0]));   // Get the sequences matrix (in packed form) 	  
	in = mxGetPr(prhs[2]);   // Get the kmer length (L)  	  
	L	= in[0]; 
	in = mxGetPr(prhs[3]);   // Get the unique flag    	  
	unique_flag = in[0]; 
	in = mxGetPr(prhs[4]);   // Get the hash flag 	  
	hash_flag = in[0]; 

	num_seqs = mxGetM(prhs[0]);    // Get number of different sequences (genes)
	in_w = (word *)(mxCalloc(num_seqs, sizeof(word))); // allocate memory to read seqs

	seqs_lens = new long[num_seqs]; // allocate the lengths array
	inn = mxGetPr(prhs[1]);   // Get the sequences original lengths in nucleotides (not packed ..) 	  	  
	//  seqs_len = (long)(inn[0]); // Give to the function the original lengths (pointer?)
	long  num_seq_lens = MAX(mxGetM(prhs[1]), mxGetN(prhs[1])); // Get number of different lengths
	for(i = 0; i < num_seqs; i++)
	{
		if(num_seq_lens>1) // read different lengths 
			seqs_lens[i] = inn[i]; // different seqs_len;    
		else   // All genes are assumed to have the same length
			seqs_lens[i] = inn[0]; // all are the same length 
		total_num_kmers += (seqs_lens[i]-L+1); // don't take from sequences ends 
		//		printf("Current seq-len=%ld, total=%ld, L=%ld num_lens=%ld\n", seqs_lens[i], total_num_kmers, L, num_seq_lens); 
		//		fflush(stdout);
	}


	//	if ( mxIsUint32(prhs[0]) != 1)
	//      mexErrMsgTxt("Input must be a 32-bit integer!!! ");


	if(nrhs > 5) // here we get also input coordinates 
	{	
		input_position_flag=1; // here we pick the positions where we extract the kmers 
		in = mxGetPr(prhs[5]);   // Get the kmer length (L)  	  
		num_coordinates = MAX(mxGetM(prhs[5]), mxGetN(prhs[5])); 
		total_num_kmers = num_coordinates; 
	}
	else
		num_coordinates = 1;


	//	printf("Obtaining %ld coordinates, %ld total kmers\n", num_coordinates, total_num_kmers); 
	//	fflush(stdout);
	kmer_seq_inds = new long[total_num_kmers];
	//	if(!input_position_flag)
	//		kmer_seq_inds = new long[total_num_kmers]; 
	kmer_position_inds = new long[num_coordinates]; 

	if(nrhs > 5) // here we get also input coordinates 
	{	
		for(i=0; i<num_coordinates; i++)
			kmer_seq_inds[i] = in[i]-1; // subtract one (matlab->c) 
		in = mxGetPr(prhs[6]);   // Get the kmer length (L)  	  
		for(i=0; i<num_coordinates; i++)
			kmer_position_inds[i] = in[i]-1; // subtract one (matlab->c) 
	}


	if(debug_prints)
	{	
		printf("Allocate seqs lengths  %ld\n", num_seqs);
		fflush(stdout);
	}




	//	printf("Inputs: kmer-length=%ld, unique-flag=%ld, hash-flag=%ld\n", L, unique_flag, hash_flag); 
	//	fflush(stdout); 



	//	printf("Allocating %ld\n", total_num_kmers);
	//	fflush(stdout);
	word **kmers = new word*[total_num_kmers]; 	//  	word **kmers = (word **) mxMalloc( total_num_kmers * sizeof(double));
	//	printf("Allocated main array kmers %ld\n", total_num_kmers);
	//	fflush(stdout);
	long *kmer_kmer_inds = new long[total_num_kmers]; 

	//	printf("Allocated main array indices %ld\n", total_num_kmers);
	//	fflush(stdout);



	//	printf("THE SEQS LEN IS : %ld\n", seqs_lens[0]);	
	//	fflush(stdout);
	in_w = (word *)(mxGetData(prhs[0]));   // Get the sequences matrix (in packed form) 	  
	seqs = new word*[num_seqs]; // read input sequences 
	for(i = 0; i < num_seqs; i++)
	{
		seqs[i] = new word[seqs_lens[i]];   // allocate memory for every sequence
		for(j = 0; j < (seqs_lens[i]-1)/half_word_size+1; j++)
			seqs[i][j] = in_w[j*num_seqs + i];  // Copy the sequences matrix (Note that they should be packed for now !!!)
	}


	// print to see that input came right
	if(debug_prints)
	{
		printf("Print input DNA: First Seq=%ld \n", seqs[0][0]); 
		for(i=0; i<num_seqs; i++)
			print_packed_dna(seqs[i], seqs_lens[i], 1);
		fflush(stdout);	
	}
	/*************************************************  Finished Read Input ***************************************************************/


	kmer_length_in_words = ceil(double(L)/double(half_word_size));
	//	printf("kmer length in words =%ld\n", kmer_length_in_words);

	for(i = 0; i < total_num_kmers; i++)
		kmers[i] = new word[kmer_length_in_words]; //   (word *) mxMalloc( kmer_length_in_words * sizeof(word));		


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
	extract_sub_kmers(L, seqs, seqs_lens, num_seqs, unique_flag, 
		kmers, kmer_kmer_inds, kmer_seq_inds, kmer_position_inds, &num_unique_kmers, 
		&num_sparse_matrix_elements, hash_flag, input_position_flag, num_coordinates);





	/*
	for(i=0; i<num_seqs; i++) 
	for(j=0; j<5; j++)
	printf("original word = %lu\n", seqs[i][j]);
	*/

	// free memory
	for(i=0; i<num_seqs; i++) 
		delete seqs[i];
	delete seqs; 
	delete seqs_lens;
	mwSize mw_total_num_kmers = total_num_kmers;
	mwSize mw_kmer_length_in_words = kmer_length_in_words;

	/*
	for(i=0; i<total_num_kmers; i++) 
	for(j=0; j<kmer_length_in_words; j++)
	printf("%ld\n", kmers[i][j]);
	*/

	//	mwSize mw_dims[2];
	//	mw_dims[0] = mw_total_num_kmers; mw_dims[1] = mw_kmer_length_in_words;


	if(debug_prints)
	{
		printf("Num. unique kmers=%ld \n", num_unique_kmers); 
		fflush(stdout); 
	}

	long return_double = 0;
	if(return_double)
	{
		plhs[0] = mxCreateDoubleMatrix(num_unique_kmers, mw_kmer_length_in_words, mxREAL); //	  plhs[0] = mxCreateNumericArray(2, mw_dims, mxUINT32_CLASS, mxREAL);
		out = mxGetPr(plhs[0]); // output the scores array 
		for(i=0; i<num_unique_kmers; i++) 
			for(j=0; j<kmer_length_in_words; j++)
				out[i+j*num_unique_kmers] = double(kmers[i][j]); // copy kmers array to output // 0*j*total_num_kmers+
	}
	else // here return integers
	{
		//		printf("Create Output Matrix[%ld,%ld] \n",num_unique_kmers, mw_kmer_length_in_words); 
		//		fflush(stdout); 
		if(word_size == 32)
			plhs[0] = mxCreateNumericMatrix(num_unique_kmers,mw_kmer_length_in_words, mxUINT32_CLASS, mxREAL);
		else
			plhs[0] = mxCreateNumericMatrix(num_unique_kmers,mw_kmer_length_in_words, mxUINT64_CLASS, mxREAL);

		//		printf("Copy to Output Matrix[%ld,%ld] \n",num_unique_kmers, mw_kmer_length_in_words); 
		//		fflush(stdout); 
		out_w = (word *)(mxGetData(plhs[0]));
		for(i=0; i<num_unique_kmers; i++) 
		{
			//			printf("Copy kmer %ld\n", i); 
			for(j=0; j<kmer_length_in_words; j++)
				out_w[i+j*num_unique_kmers] = (kmers[i][j]); // copy kmers array to output // 0*j*total_num_kmers+
		}
	}


	for(i = 0; i < total_num_kmers; i++)
		delete kmers[i]; 			  
	delete kmers; //   mxFree(kmers);

	//  printf("one Shift 32 is: %lu, and %lu\n", 1UL << 32, 143241231 + (1UL << 32)); 
	plhs[1] = mxCreateDoubleMatrix(num_sparse_matrix_elements, 2, mxREAL); //		plhs[1] = mxCreateNumericMatrix(total_num_kmers, 1, mxINT32_CLASS, mxREAL);

	out = mxGetPr(plhs[1]); // output the scores array 
	for(i=0; i<num_sparse_matrix_elements; i++) 
	{
		//		for(j=0; j<2; j++)
		out[i+0] = double(kmer_kmer_inds[i]+1); // copy scores array to output. add 1 (indices in matlab start at 1) 
		out[i+num_sparse_matrix_elements] = double(kmer_seq_inds[i]+1); // copy scores array to output. add 1 (indices in matlab start at 1) 
	}
	//	for(i = 0; i < total_num_kmers; i++)
	//		delete kmer_inds[i]; 
	delete kmer_kmer_inds; 
	delete kmer_seq_inds; 
	delete kmer_position_inds;

	//	for(i=0; 0<1; i++) // infinite loop 
	//		j=i; // sleep(1000); 



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
//  function : extract_sub_kmers
// 
//  input    : L - length of kmers to extract (in nucleotides)
//			   seqs - array of sequences (packed form 2-dimensional array)
//			   seqs_lens - array of sequences lengths (num. of nucleotides in each sequence) 
//			   hash_flag - do we use hash in program? (currently always use hash)						   
//	
//
//  output   : kmers - array of unique kmers (in packed form 2-dim array) 
//			   kmer_kmer_inds - array of indices for each kmer. Index in the list of unique kmers
//			   kmer_seq_inds - array of indices for each kmer. Index is in the list of sequences (species) 
//			   kmer_position_inds - array of indices for each kmer. Index is the position in the sequences (in each species)
//			   num_unique_kmers - total number of unique kmers
//			   num_sparse_matrix_elements - total number of (non-unique) kmers to be used in sparse matrix
//
//  purpose  : Scan sequences and extract all subsequences of length L. 
//			   Use a hash table to speed the process.   
//			   
//
////////////////////////////////////////////

long extract_sub_kmers(long L, word *seqs[MAX_NUM_SEQS], long *seqs_lens, long num_seqs, long unique_flag, 
					   word *kmers[MAX_L], long *kmer_kmer_inds, long *kmer_seq_inds, long *kmer_position_inds, 					   
					   long *num_unique_kmers, long *num_sparse_matrix_elements, long hash_flag, 
					   long input_position_flag, long num_coordinates)
{
	long i, j, k;
	long i_seq, j_position; 
	long num_kmers_enum;
	long half_word_size = 4*sizeof(word);
	long word_size = 2*half_word_size;

	long num_words_in_kmer = ceil(double(L) / double(half_word_size)); 
	long num_bytes_in_kmer = ceil(double(L) / 4.0); // temp - a bit wasteful but easier to make sure it's correct 
	long num_bytes_in_kmer_ceil = 4 * ceil(double(L) / 16.0); // temp - a bit wasteful but easier to make sure it's correct 

	long last_word_in_kmer = (L-1)%half_word_size+1; 
	long kmer_ctr=0; // count which unique kmer are we at 
	long index_ctr=0; // count which index are we at the sequence (non-unique kmer) 
	struct SeqData *found_hash_ptr;
	struct SeqData *temp_ptr;
	long current_idx_in_hash;
	word current_kmer[MAX_L]; // save one kmer 
	long debug_prints = 0; 
	long rand_int2;
	long add_int = 0; 
	word ones_mask;
	if(word_size == 32)
		ones_mask = 0xFFFFFFFF;
	else // assume it's 64 
		ones_mask = 0xFFFFFFFFFFFFFFFF;

	if(debug_prints)
	{
		printf("ONES_MASK=%lx\n", ones_mask); 
		printf("ONES_MASK_DEC=%lu\n", ones_mask); 
		printf("SIZEOFLONG=%ld (bits)\n", 8*sizeof(long)); 
		printf("kmer_len=%ld, Num. Words in kmer=%ld, num. bytes=%ld\n", L, num_words_in_kmer, num_bytes_in_kmer); 
	}

	long num_seqs_enum = num_seqs;
	//	printf("Input positions = %ld\n", input_position_flag); 
	if(input_position_flag)
	{
		num_seqs_enum = num_coordinates;
		if(debug_prints)
		{
			for(i=0; i<num_coordinates; i++)
				printf("Pos %ld=%ld\n", i, kmer_position_inds[i]); 

			printf("Run on %ld Sequences\n", num_seqs_enum); 
			for(i = 0; i < num_seqs_enum; i++) // loop on all sequences 
			{
				printf("Before loop position ind %ld= %ld\n", i, kmer_position_inds[i]);
				fflush(stdout); 
			}
		}

	}


	for(i = 0; i < num_seqs_enum; i++) // loop on all sequences 
	{

		if(debug_prints)
		{
			printf("Run on seq %ld\n", i); 
			fflush(stdout); 
		}
		if(input_position_flag)
		{
			i_seq = kmer_seq_inds[i]; 
			num_kmers_enum = 1;
		}
		else
		{
			i_seq = i;
			num_kmers_enum = seqs_lens[i]-L+1;
		}
		if(debug_prints)
		{		printf("Cur seq ind: %ld\n", i_seq);
		fflush(stdout); 
		}
		for(j = 0; j < num_kmers_enum; j++) // loop on all possible kmers in a sequence 
		{
			if(input_position_flag)
			{
				j_position = kmer_position_inds[i]; 
				//				printf("Cur position ind: %ld\n", j_position);
				//				fflush(stdout); 
			}
			else
				j_position = j;

			//			printf("Use words: %lu %lu, kmer_ctr=%ld \n", seqs[i][j/half_word_size], seqs[i][j/half_word_size+1], kmer_ctr);
			for(k = 0; k < num_words_in_kmer; k++) // full words (get 16 nucleotides each time)
			{
				//				printf("run word %ld, i_seq=%ld, j_position=%ld, ind_read=(%ld,%ld,%ld)", k, i_seq, 
				//					j_position, j_position/half_word_size+k, j_position/half_word_size+1+k); 
				//				fflush(stdout);
				current_kmer[k] = ((seqs[i_seq][j_position/half_word_size+k] >> (2*(j_position%half_word_size)))&ones_mask) + 
					(((seqs[i_seq][j_position/half_word_size+1+k]&((1UL << (2*(j_position%half_word_size)))-1)) << 
					(2*(half_word_size-(j_position%half_word_size))))&ones_mask); // merge two half words from both adjacent words 
				//				printf(" word_kmer=%ld, ", current_kmer[k]); 
				//				fflush(stdout);


				//				printf("LAST WORD LEN = %ld\n", last_word_in_kmer); 
				if((k == num_words_in_kmer-1) && (last_word_in_kmer < half_word_size)) 
				{
					current_kmer[k] = current_kmer[k]&((1UL << (2*last_word_in_kmer))-1); // take only L nucleotides from last word 
					//				printf(" Last Word!!!\n"); 
					//				fflush(stdout);
				}

				//				kmers[kmer_ctr][k] = kmers[kmer_ctr][k]&((1UL << (2*last_word_in_kmer))-1); // take only L nucleotides
			}			

			if(debug_prints)
			{
				printf("Extracted Kmer: (L=%ld)\n", L); 
				print_packed_dna(current_kmer, L, 1); // print extracted kmer 
				fflush(stdout); 
			}

			if( (!unique_flag) || (!hash_flag) ) // here we just concatenate 
			{							
				if(debug_prints)
				{
					printf("Problem! Not Unique!\n");
					fflush(stdout); 

					printf("seq kmer[%ld]=%ld\n", kmer_ctr, i_seq);
					fflush(stdout);
					printf("index_ctr=%ld\n", index_ctr); 
					fflush(stdout);
					printf("num_bytes_to_copy=%ld\n", num_bytes_in_kmer_ceil); 
					fflush(stdout);
				}
				kmer_kmer_inds[index_ctr] =  kmer_ctr; //  index where kmer is in the list of all kmers 
				kmer_seq_inds[index_ctr++] = i_seq; // index of species for this kmer
				memcpy(kmers[kmer_ctr++], current_kmer, num_bytes_in_kmer_ceil); // copy kmer into output list 
				//				printf("EEExtracted Kmer: (L=%ld)\n", L); 
				//					print_packed_dna(current_kmer, L, 1); // print extracted kmer 
				//				fflush(stdout); 
				/**
				if(problematic_kmer)					// print kmer again
				{
				printf("Put prob. kmer in list place=%ld: \n", kmer_ctr); 						
				print_packed_dna(kmers[kmer_ctr-1], 32, 1);

				}
				**/

				//				kmer_seq_inds[kmer_ctr++] = i_seq; // record which sequence did the kmer come from 	(there's only one index here)		kmer_ctr++; 
			}
			else // keep a unique list by using Hash 
			{
				//					if(debug_prints)
				//					{
				//						printf("Find in Hash!!!\n"); 
				//						fflush(stdout); 
				//					}
				if(add_int)
				{
					rand_int2 = rand(); 
					HASH_FIND_INT(hTable, &rand_int2, found_hash_ptr);
				}
				else
					HASH_FIND(hh, hTable, current_kmer, num_bytes_in_kmer, found_hash_ptr);  

				//					found_hash_ptr = NULL; 

				//					if(debug_prints)
				//					{
				//						printf("Find in Hash22222\n"); 
				//						fflush(stdout); 
				//					}
				if(found_hash_ptr == NULL) // found a new kmer 
				{
					if(i%50 == 0)
					{
						if(debug_prints)
						{
							printf("Not found in Hash New ID=%ld, seq=%ld pos=%ld, kmer=%ld\n", kmer_ctr, i, j, current_kmer[0]); 
							fflush(stdout); 
						}
					}

					temp_ptr = (struct SeqData *) calloc(1,sizeof(struct SeqData));	
					memcpy(kmers[kmer_ctr], current_kmer, num_bytes_in_kmer); // copy kmer into output list 
					memcpy(temp_ptr->seq, current_kmer, num_bytes_in_kmer); // copy kmer into hash table 
					temp_ptr->seqpos = kmer_ctr; // copy index (what is this?)

					if(add_int)
					{
						temp_ptr->rand_int = rand(); 
						HASH_ADD_INT(hTable, rand_int, temp_ptr); // add value to Hash table 
					}
					else
					{
						HASH_ADD(hh, hTable, seq, num_bytes_in_kmer, temp_ptr); // add value to Hash table 
					}
					kmer_kmer_inds[index_ctr] = kmer_ctr++; // save index of this kmer in the kmers long list
				}
				else // kmer already appears in hash table. Just update indices 
				{	
					if(i%50 == 0)
					{
						if(debug_prints)
						{
							printf("Found in Hash Before, seq=%ld pos=%ld, kmer=%ld\n", i, j, current_kmer[0]); 
							fflush(stdout); 
							printf("Found in Hash OLD ID=%ld\n", found_hash_ptr->seqpos); 
							fflush(stdout); 
						}
					}
					kmer_kmer_inds[index_ctr] = found_hash_ptr->seqpos; // index where kmer was found in hash.  save index of this kmer in the kmers long list 
					//						if(debug_prints)
					//						{
					//							printf("Found Hash After\n");
					//							fflush(stdout); 
					//						}
				}
				kmer_seq_inds[index_ctr++] = i_seq; // save index of species for this kmer 
			} // end if (non unique or non hash) 

			//			printf("Ind: %ld\n", kmer_inds[kmer_ctr-1]); 
			//			kmers[kmer_ctr][num_words_in_kmer-1] = // last word may be part 
			//				seqs[i][j/HALF_WORD_SIZE+1+k]
		}
	}
	/**/
	if( (!hash_flag) && unique_flag ) // perform sort and unique 
	{

//		printf("No hashing used \n"); 
		//		word *tmp_kmers;
		//		tmp_kmers = new word[kmer_ctr]; 
		long *tmp_kmer_perm_indices;
		tmp_kmer_perm_indices = new long[kmer_ctr]; 
		long *kmer_keep_inds; kmer_keep_inds = new long[kmer_ctr]; 


		if(debug_prints)
		{
			printf("Before QSort:\n");
			for(i=0; i<10/*kmer_ctr*/; i++)
			{
				printf("\n i=%ld kmer= ", i); 		
				for(j=0; j<num_words_in_kmer; j++)
					printf(" %lu, ", kmers[i][j]); 		
			}
		}
		/**
		for(i=0; i<kmer_ctr; i++)
		if(kmers[i][0] == 460002223925673672)
		{
		printf("Before Qsort Found KMER!! pos=%ld\n", i); 
		}
		**/
		DoQuicksort(kmers, kmer_ctr, num_words_in_kmer, num_bytes_in_kmer, tmp_kmer_perm_indices); // sort the ORIGINAL KMERS 		
		/**		
		for(i=0; i<kmer_ctr; i++)
		if(kmers[i][0] == 460002223925673672)
		{
		printf("Before Qsort Found KMER!! pos=%ld\n", i); 
		}

		**/
		if(debug_prints)
		{
			printf("\n\nAfter QSort:\n");
			for(i=0; i<kmer_ctr; i++)
			{
				//				printf("\n i=%ld kmer= ", i); 		
				//				for(j=0; j<num_words_in_kmer; j++)
				//					printf(" %lu, ", kmers[i][j]); 		
			}
			printf("LEN WORD USED=%ld\n", num_words_in_kmer);
		}
		//		for(i=0; i<kmer_ctr; i++)
		//			printf("%ld\n", kmers[i][0]); 		


		DoOrder(kmer_seq_inds, kmer_ctr, tmp_kmer_perm_indices);
		delete tmp_kmer_perm_indices;
		index_ctr = kmer_ctr; 
		kmer_ctr = -1; // first time this is incremented to zero  
		for(i=0; i<index_ctr; i++) // now run over all kmers and do unique 
		{
			//			printf("i=%ld kmer=%ld. ", i, tmp_kmers[i]); 
			if((i>0) && (!memcmp(kmers[i], kmers[i-1], num_bytes_in_kmer)))	 // same as previous 
			{
				kmer_kmer_inds[i] = kmer_ctr; // index of this kmer in list of all kmers 
				//				kmer_seq_inds[i] = tmp_kmer_inds[1][i]; // index of species for this kmer 
				if(debug_prints)
				{
					//					printf("Same kmers again!!!!\n"); 
					//					fflush(stdout); 
				}

			}
			else // found a new value 
			{
				kmer_ctr++;
				kmer_keep_inds[kmer_ctr] = i; 
				//				memcpy(kmers[kmer_ctr], kmers[i], num_bytes_in_kmer); 
				//				kmers[kmer_ctr][0] = tmp_kmers[i]; // copy new value 
				kmer_kmer_inds[i] = kmer_ctr; // index of this kmer in list of all kmers 
				//				kmer_inds[i][1] = tmp_kmer_inds[1][i]; // index of species for this kmer 
				if(debug_prints)
				{
					//					printf("Incrementing kmers!!!!\n"); 
					//					fflush(stdout); 
				}
			}
		} // finish loop on all kmers
		kmer_ctr++; // true number of found kmers 
		//		printf("NUM UNIQUE KMERS IN EXTRACT FUNCTION: %ld\n", kmer_ctr); 
		for(i=0; i<kmer_ctr; i++) // copy again into same array
		{
			if(kmer_keep_inds[i] > i)
				memcpy(kmers[i], kmers[kmer_keep_inds[i]], num_bytes_in_kmer);
		}

		delete kmer_keep_inds;
		//		DoUnique(kmers); 
	}
	else
	{
		if(debug_prints)
		{
		if(hash_flag)
			printf("Used hash \n"); 
		else
			printf("No unique needed \n"); 
		}
	}
	/**/	




	/*
	for(i=0; i<kmer_ctr; i++)
	printf("kmer[%ld]=%lu\n", i, kmers[i][0]); 
	printf("Found %ld unique kmers\n", kmer_ctr);  
	*/

	num_unique_kmers[0] = kmer_ctr; 
	num_sparse_matrix_elements[0] = index_ctr;

	delete_all(); // free hash memory 
	if(hash_flag)
		HASH_CLEAR(hh,hTable); // free memory 

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
	word ones_mask=0xF;
	if(half_word_size == 32)
		ones_mask = 0x1F; 


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
