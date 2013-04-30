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
/// New Functions For Yuval's Project - 11.11.2003
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])
{

	int nRows, nCols, i, j, K;

	double weights[4][MAX_L];
	float weights_single[4][MAX_L];	
	long L; 
	word **seqs;
	double **scores;
	float **scores_single;
	long *seqs_lens;   // New : There may be different sequences lengths
	long num_seqs; 
	long seqs_len; 
	double threshold;
	long half_word_size = 4*sizeof(word); // half word size in bits (=word size in nucleotides)
	long word_size = 2*half_word_size;

	double *in, *inn, *out;
	float *in_single, *out_single;	

	/* Check for proper number of arguments. Should be five 5 (last one is double flag) */
	/* Vars are: pwms, seqs, seqs_lens, -10, is_pwm_double, is_seqs_double */
	if( (nrhs < 5) || (nrhs > 6) ) {
		mexErrMsgTxt("Usage: [pwms] [seqs] [seqs_lens] [-10] [is_pwm_double] (optional) [is_seqs_double]");
	}

	in = mxGetPr(prhs[4]);   // Get the double flag matrix
	long is_pwm_double = long(in[0]);
	long is_seqs_double;
	if(nrhs == 6) // here we set the input type (word or double)
	{
		in = mxGetPr(prhs[5]);   // Get the double flag matrix
		is_seqs_double = long(in[0]);
	}
	else 
		is_seqs_double = 0; 

	if(is_pwm_double)
	{
		in = mxGetPr(prhs[0]);   // Get the weight matrix
		L = mxGetN(prhs[0]);     // Get the weight matrix's length  
		for(j = 0; j < L; j++)
			for(i = 0; i < 4; i++)
				weights[i][j] = in[j*4+i];   // Copy the weight matrix (pwm)
	}
	else
	{
		in_single = (float *)(mxGetData(prhs[0]));		// Get the weight matrix (singles)
		L = mxGetN(prhs[0]);     // Get the weight matrix's length  
		for(j = 0; j < L; j++)
			for(i = 0; i < 4; i++)
				weights_single[i][j] = in_single[j*4+i];   // Copy the weight matrix (pwm)
	}

	/**
	printf("At the begining weights : \n");
	for(i = 0; i < 4; i++)
	{
	printf("\n");
	for(j = 0; j < L; j++)
	printf("%lf   ", weights[i][j]);
	}
	fflush(stdout);
	/**/

	in = mxGetPr(prhs[1]);   // Get the sequences matrix	  	  
	num_seqs = mxGetM(prhs[1]);    // Get the number of sequences (genes)
	word *in_w = (word *)(mxGetData(prhs[0]));   // Get the sequences matrix (in packed form) 	  


	inn = mxGetPr(prhs[2]);   // Get the sequences original length in nucleotides (not packed ..) 	  	  
	seqs_len = long(inn[0]); // Give to the function the original lengths
	//  printf("THE SEQS LEN IS : %ld\n", seqs_len);	

	seqs = new word*[num_seqs];
	for(i = 0; i < num_seqs; i++)
	{
		seqs[i] = new word[long(1+(seqs_len+1)/half_word_size)];   // allocate memory for every sequence (why so much??) 
		if(is_seqs_double) // new!!! 2 means that 
		{
			for(j = 0; j < (seqs_len-1)/half_word_size+1; j++)
				seqs[i][j] = in[j *  num_seqs + i  /*i*num_seqs+j*/];  // Copy the sequences matrix (Note that they should be packed for now !!!)
		}
		else
		{
			for(j = 0; j < (seqs_len-1)/half_word_size+1; j++)
				seqs[i][j] = in_w[j *  num_seqs + i  /*i*num_seqs+j*/];  // Copy the sequences matrix (Note that they should be packed for now !!!)
		}
	}

	in = mxGetPr(prhs[3]);   // Get the  threshold	  	  
	threshold	= in[0]; // currently meaningless

	if(is_pwm_double)
	{
		scores = new double*[num_seqs];
		for(i = 0; i < num_seqs; i++)
			scores[i] = new double[seqs_len-L+1];  // allocate memory for the scores
	}
	else
	{
		scores_single = new float*[num_seqs];
		for(i = 0; i < num_seqs; i++)
			scores_single[i] = new float[seqs_len-L+1];  // allocate memory for the scores
	}

	seqs_lens = new long[num_seqs]; // allocate the lengths array
	for(i = 0; i < num_seqs; i++)
		seqs_lens[i] = seqs_len;    // All genes are assumed to have the same length

	///  printf("Paramters For calc_all_sites_scores : L = %ld num_seqs = %ld seqs_lens[0] = %ld thresh = %ld\n", 
	///	  L, num_seqs, seqs_lens[0], threshold);
	///  fflush(stdout);

	if(is_pwm_double)
	{
		calc_all_sites_scores(weights, L, seqs, 
			scores, num_seqs, seqs_lens);
		plhs[0] = mxCreateDoubleMatrix(num_seqs, seqs_len-L+1, mxREAL);
		out = mxGetPr(plhs[0]); // output the scores array 
		for(i=0; i<num_seqs; i++) 
			for (j=0; j<seqs_len-L+1; j++)
				out[j*num_seqs+i] = scores[i][j]; // copy scores array to output
	}
	else
	{
		calc_all_sites_scores_single(weights_single, L, seqs, 
			scores_single, num_seqs, seqs_lens);
		plhs[0] = mxCreateNumericMatrix(num_seqs, seqs_len-L+1, mxSINGLE_CLASS, mxREAL);
		out_single = (float *)(mxGetData(plhs[0])); // output the scores array 
		for(i=0; i<num_seqs; i++) 
			for (j=0; j<seqs_len-L+1; j++)
				out_single[j*num_seqs+i] = scores_single[i][j]; // copy scores array to output
	}



	// free memory
	for(i=0; i<num_seqs; i++) 
		delete seqs[i];

	if(is_pwm_double)
	{
		for(i=0; i<num_seqs; i++) 
			delete scores[i];
		delete scores;
	}
	else
	{
		for(i=0; i<num_seqs; i++) 
			delete scores_single[i];
		delete scores_single;
	}

	delete seqs; 
	delete seqs_lens;


}
//#endif // MEX_TO_MATLAB
