#include <stdlib.h>
#include <stdio.h>
//#include <iostream.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "general.h"
#include "dna_utils_new.h"
#include "markov.h"



///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For Yuval's Project - 04.07.2005
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
  int nRows, nCols, i, j, K;

  double weights[4][MAX_L];
  long L; 
  word **seqs;
  double **scores;
  long num_seqs; 

  double threshold;
  
  double *in, *inn, *out;
  long *seqs_lens; // If we have many sequences with different lengthes   

//	printf("STARTED MEX FUNN\n");
//	fflush(stdout);

  /* Check for proper number of arguments. Should be Seven 7 */
  if(nrhs != 4) {
    mexErrMsgTxt("Usage: [dat] [K]");
  }
  
  in = mxGetPr(prhs[0]);   // Get the weight matrix
  L = mxGetN(prhs[0]);     // Get the weight matrix's length
  
  for(j = 0; j < L; j++)
	 for(i = 0; i < 4; i++)
		  weights[i][j] = in[j*4+i];   // Copy the weight matrix

// Note !!! here the sequences come in one vector, i.e. matrix of size one !!! 
  in = mxGetPr(prhs[1]);   // Get the sequences matrix	(actually a vector !!)  	  

  inn = mxGetPr(prhs[2]);   // Get the sequences original lengths (not packed ..) 	  	  
  num_seqs = mxGetN(prhs[2]);    // Get the number of sequences (genes)

  seqs_lens = new long[num_seqs]; // New ! allocate the lengths array
  for(i = 0; i < num_seqs; i++)
	  seqs_lens[i] = long(inn[i]);    // Get the different lengths of the genes. 

//  printf("There are %ld Seqs in calc_scores ! THE SEQS LEN IS :\n", num_seqs);	
//  for(i = 0; i < num_seqs; i++)
//	  printf("[%ld] --> %ld <--\n", i, seqs_lens[i]); 


  seqs = new word*[num_seqs];
  long ind = 0; // index of the locations in the sequences
  for(i = 0; i < num_seqs; i++)
  {
	  seqs[i] = new word[MAX((seqs_lens[i]-1)/16+1,0)]; //// MAX(seqs_lens[i]-1)/16+1,0)  // allocate memory for every sequence
	  for(j = 0; j < (seqs_lens[i]-1)/16+1; j++)
			seqs[i][j] = in[ind + j ];  // Copy the sequences matrix (Note that they should be packed for now !!!)
	  ind += MAX(((seqs_lens[i]-1)/16+1),0); //// MAX(((seqs_lens[i]-1)/16+1),0);
  }

  // Create a matrix with the maximal seqs_len
  long sum_seqs_len = 0; 
  for(i=0; i<num_seqs; i++) 
	  sum_seqs_len += (MAX(seqs_lens[i]-L+1,0)); //// MAX((seqs_lens[i]-L+1), 0);


  // Set the memory for the scores : 
  scores = new double*[num_seqs];
  for(i = 0; i < num_seqs; i++)
	scores[i] = new double[seqs_lens[i]-L+1]; //// double[MAX(seqs_lens[i]-L+1,0)];  // allocate memory for the scores


  in = mxGetPr(prhs[3]);   // Get the  threshold	  	  
  threshold	= in[0]; // currently meaningless


  calc_all_sites_scores(weights, L, seqs, 
  						scores, num_seqs, seqs_lens);

  // Here's the problem : We need to give the output in one array and not a matrix !!!!
  plhs[0] = mxCreateDoubleMatrix(1, sum_seqs_len, mxREAL);
  out = mxGetPr(plhs[0]); // output the scores array 
  ind=0;
  for(i=0; i<num_seqs; i++) 
  {
	for (j=0; j<(seqs_lens[i]-L+1); j++)  // Here we go up to the original length !!!! 	
		out[ind+j] = scores[i][j]; // copy scores array to output
	ind += (MAX(seqs_lens[i]-L+1,0)); ////	ind += MAX(seqs_lens[i]-L+1,0);  // new : avoid negative increments !!!! 	
  }

  // free memory
  for(i=0; i<num_seqs; i++) 
  {
	  delete scores[i];
	  delete seqs[i];
  }
  delete scores;
  delete seqs; 
  delete seqs_lens;

}
//#endif // MEX_TO_MATLAB


// Print An array In Binary
void print_bin(word *a, long L)
{
	long i;

	for(i = 0; i < L; i++)
		printf("%d", (a[i/WORD_SIZE] >> (i%WORD_SIZE))&1);
	printf("\n");

}



// Print A number In 2-bit digits
void print_quad(word *a, long L)
{
	long i;

	for(i = 0; i < L; i++)
		printf("%d", (a[i/HALF_WORD_SIZE] >> (2*(i%HALF_WORD_SIZE)))&3);
	printf("\n");

} 



// Print A number In 2-bit RNA digits
void print_rna(word *a, long L)
{
	long i;
	word digit;

	for(i = 0; i < L; i++)
	{
		digit = (a[i/HALF_WORD_SIZE] >> (2*(i%HALF_WORD_SIZE)))&3;

		if(digit == 0)
			printf("a");
		if(digit == 1)
			printf("c");
		if(digit == 2)
			printf("g");
		if(digit == 3)
			printf("u");

	}
	printf("\n");

} 

