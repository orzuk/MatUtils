#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "general.h"
#include "hmm_chrom_funcs.h"
//#include "dna_utils_new.h"
//#include "markov.h"







///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For HMM Chromosome project - 11.07.2004
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB


// The function prepares all the data needed for the function TrainModelEM.
// The structures for the hmm and the data are build. Other parameters
// such as iterations, tolerance etc. are also supplied
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
	int nRows, nCols, i, j, K, t, max_iters, num_starts;

	long seq_len, L; 
	long x_dim, y_dim, place_flag, miss_data, fold_change_flag;
	double tolerance;
  
	double *in, *out;

	double best_model_score;


	hmm_model hmm;   // a struct containing the hmm parameters

	hmm_data data;  // a struct containing all the data needed


	 /* Check for proper number of arguments. Should be Fourteen 14 */
	if(nrhs != 14) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. Expression vector
	2. Chromosomal location vector
	3. Missing data flag (says if to use locations ... )
	4. X dimension required
	5. Y dimension (mixtures) required
	6. place flag (if to use different distribution for each place)
	7. fold change flag (if to bound the mus by one) 
	8. mean vec (currently it is given and not learned)
	9. std vec (currently it is given and not learned)
	10. UpperBounds matrix (bounds on entries of the transition matrix)
	11. bound flag (says if to use upper bounds)
	12. number of EM iterations
	13. Number of EM starting points 
	14. EM tolerance 
	*/


	in = mxGetPr(prhs[0]);   // Get the expression vector
	seq_len = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the expression vector length

	data.y_vec = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy expression vector
		data.y_vec[t] = in[t];  
	

	in = mxGetPr(prhs[1]);   // Get the location vector
	L  = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));    // Get the number of sequences (genes)
	if(L != seq_len)
		mexErrMsgTxt("Usage: Exp. and Loc. vector length must be the same");

	data.loc_vec = new long[seq_len];
	for(t = 0; t < seq_len; t++)   // copy location vector
		data.loc_vec[t] = long(in[t]/NUM_BASEPAIRS_PER_HMM_TRANS);  

	data.loc_diff_vec = new long[seq_len];
	data.loc_diff_vec[0] = 0; 
	for(t = 1; t < seq_len; t++)   // copy location vector
		data.loc_diff_vec[t] = MIN(MAX(data.loc_vec[t]-data.loc_vec[t-1]-1, 0), MAX_POWER_OF_TRANS_MATRIX-2); 
	
/*  // Print the distances 
	for(t = 0; t < seq_len; t++) 
	{
		if(data.loc_diff_vec[t] != 0)
			printf("%ld -> %ld \n", t, data.loc_diff_vec[t]);
	}
	printf("\n");
*/	

	in = mxGetPr(prhs[2]);   // Get the  missing data flag
	miss_data = long(in[0]); 

	in = mxGetPr(prhs[3]);   // Get the x-dimension 
	x_dim = long(in[0]);     // Copy the x-dimension

	in = mxGetPr(prhs[4]);   // Get the y-dimension 
	y_dim = long(in[0]);     // Copy the y-dimension
	

	// here insert the place 
	in = mxGetPr(prhs[5]);   // Get the place flag 
	place_flag = long(in[0]);     // Copy the place flag

	// here insert the fold change 
	in = mxGetPr(prhs[6]);				// Get the fold-change flag 
	fold_change_flag = long(in[0]);     // Copy the fold-change flag
	


	if(place_flag)
	{
		
		// allocate memory 
		for(i = 0; i < x_dim; i++)		
			for(j = 0; j < y_dim; j++)
			{
				hmm.place_gaussian_mu[i][j] = new double[seq_len];
				hmm.place_gaussian_sigma[i][j] = new double[seq_len];
			}


	/***/

		// Now copy
		in = mxGetPr(prhs[7]);   // Get the place means
		long arr_size = mxGetM(prhs[7]) * mxGetN(prhs[7]);
		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
			{
				for(t = 0; t < seq_len; t++)
				{
					if(t+seq_len*j+seq_len*y_dim*i >= arr_size)
					{
						printf("ERROR !! TOO LARGE INDEX %ld\n", t+seq_len*j+seq_len*y_dim*i);
						break;
					}
					hmm.place_gaussian_mu[i][j][t] = in[t+seq_len*j+seq_len*y_dim*i]; // Copy the place means
				}
			}

		printf("Copied mu\n");

		// Now copy
		in = mxGetPr(prhs[8]);   // Get the place sigmas
		arr_size = mxGetM(prhs[8]) * mxGetN(prhs[8]);

		for(i = 0; i < x_dim; i++)	
			for(j = 0; j < y_dim; j++)
			{
//				printf("DO SIGMA (I J ) %ld %ld \n", i, j);
				for(t = 0; t < seq_len; t++)
				{	
					if(t+seq_len*j+seq_len*y_dim*i >= arr_size)
					{
						printf("ERROR !! TOO LARGE INDEX %ld\n", t+seq_len*j+seq_len*y_dim*i);
						break;
					}
					hmm.place_gaussian_sigma[i][j][t] = in[t+seq_len*j+seq_len*y_dim*i]; // Copy the place sigmas
				}
			}

		printf("Copied sigma\n");

 /***/
	}
	/////////////////
		


	// Now determine the upperbounds
	in = mxGetPr(prhs[9]);   // Get the upperbounds
	for(i = 0; i < x_dim; i++)	
		for(j = 0; j < x_dim; j++)
			hmm.M_upperbounds[i][j] = in[i*x_dim+j];


	in = mxGetPr(prhs[10]);   // Get the upperbounds flag
	hmm.use_bounds = long(in[0]);


	in = mxGetPr(prhs[11]);   // Get the iterations 
	max_iters = long(in[0]); // Copy the iteations

	in = mxGetPr(prhs[12]);   // Get the number of starting points 
	num_starts = long(in[0]); // Copy the number of starting points

	in = mxGetPr(prhs[13]);   // Get the tolerance for stopping iterations 
	tolerance = in[0];		 // Copy the tolerance for stopping iterations
 

	data.seq_len = seq_len; // this does in the data
	data.y_type = CONTINUOUS;
	data.miss_data = miss_data;
	data.dont_miss_prob = 0.5; // not important (only for simulations)


	hmm.miss_data = miss_data; // double .. 
	hmm.x_dim = x_dim;   
 	hmm.y_dim = y_dim;
	hmm.y_type = CONTINUOUS; // Gaussian !!!!  Should add a flag saying if we want discrete outputs as well !
	hmm.place_flag = place_flag; /// place_flag; // currently no placing ? 
	hmm.fold_change_flag = fold_change_flag;   // copy the fold-change flag 
	
/**
	printf("Sending to Function : \n");
	printf("len %ld iters %ld starts %ld tol %lf \n", seq_len, max_iters, num_starts, tolerance);

	for(t = 0; t < seq_len; t++)
		printf("%lf ,", y_vec[t]);
	printf("\n");

	PrintDoubleVec(y_vec, seq_len);
	printf("\n---\n\n");
/**/
  
//	printf("Calling EM !!! USE BOUNDS : %ld \n", hmm.use_bounds);

	// Call c function to do the job for you ..
	TrainModelEM(&hmm, &data, max_iters, num_starts, tolerance, &best_model_score);



//	printf("After Training the trained model is : \n");
	PrintModel(&hmm);


	// Now copy output from the HMM structure  : 
	plhs[0] = mxCreateDoubleMatrix(1, x_dim, mxREAL);  // The initial vector PI
	out = mxGetPr(plhs[0]); // output PI
	for(i=0; i<x_dim; i++) 
		out[i] = hmm.PI[i];

	plhs[1] = mxCreateDoubleMatrix(x_dim, x_dim, mxREAL);  // The transition matrix M
	out = mxGetPr(plhs[1]); // output M
	for(i=0; i<x_dim; i++) 
		for(j=0; j<x_dim; j++)
			out[j*x_dim+i] = hmm.M[i][j];

	plhs[2] = mxCreateDoubleMatrix(x_dim, y_dim, mxREAL);  // The mixture matrix N
	out = mxGetPr(plhs[2]); // output N
	for(i = 0; i < x_dim; i++) 
		for(j = 0; j < y_dim; j++)
			out[j*x_dim+i] = hmm.N[i][j];

	plhs[3] = mxCreateDoubleMatrix(x_dim, y_dim, mxREAL);  // The mean matrix mu
	out = mxGetPr(plhs[3]); // output mu
	for(i = 0; i < x_dim; i++) 
		for(j = 0; j < y_dim; j++)
			out[j*x_dim+i] = hmm.MU[i][j];

	plhs[4] = mxCreateDoubleMatrix(x_dim, y_dim, mxREAL);  // The mixture matrix SIGMA
	out = mxGetPr(plhs[4]); // output SIGMA
	for(i = 0; i < x_dim; i++) 
		for(j = 0; j < y_dim; j++)
			out[j*x_dim+i] = hmm.SIGMA[i][j];


	plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);  // The log-score
	out = mxGetPr(plhs[5]); // output SIGMA
	out[0] = best_model_score;




/***
	plhs[0] = mxCreateDoubleMatrix(num_seqs, seqs_len-L+1, mxREAL);
	out = mxGetPr(plhs[0]); // output the scores array 
	for(i=0; i<num_seqs; i++) 
		for (j=0; j<seqs_len-L+1; j++)
			out[j*num_seqs+i] = scores[i][j]; // copy scores array to output

***/		

  // free memory
	delete data.y_vec;
	delete data.loc_vec; 
	delete data.loc_diff_vec;

	if(place_flag)
	{	

		// delete memory 
		for(i = 0; i < hmm.x_dim; i++)		
			for(j = 0; j < hmm.y_dim; j++)
			{
				delete hmm.place_gaussian_mu[i][j];
				delete hmm.place_gaussian_sigma[i][j];
			}
	}


}
//#endif // MEX_TO_MATLAB




