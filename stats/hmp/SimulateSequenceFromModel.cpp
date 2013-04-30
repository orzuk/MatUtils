#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "general.h"
#include "hmm_chrom_funcs.h"


long PrintModel2(hmm_model *hmm);

///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For HMM Chromosome project - 11.07.2004
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
	int nRows, nCols, i, j, K, t, max_iters, num_starts;

	long seq_len, L; 
	long x_dim, y_dim; // , place_flag, miss_data;
	double tolerance;
  
	double *in, *out;

//	long *y_vec_int; // dummy !!! 
//	double *y_vec; 
//	double *loc_vec;  // location on chromosome vector - Currently not used !!!! 
	
	hmm_model hmm;  // containts the model

	hmm_data data; // contains the data

	double *b_cond_tabs[MAX_X_VALS]; // helpful conditional tables
	double *alpha[MAX_X_VALS];
	double *beta[MAX_X_VALS];
	double *gamma[MAX_X_VALS];
	double *gammagamma[MAX_X_VALS][MAX_Y_VALS];
	double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS];


	long *opt_x_vec_out;  // optimal viterbi path

	double *scale, *scale_cum, *scale_exp;
	double cur_score;


	
	long miss_data = 0; // Don't miss data for now !!!! 
	long place_flag = 0; // No place flad for now !
	long use_bounds = 0; // no bounds on hmm probabilities for now


	
	/* Check for proper number of arguments. Should be Seven 7 */
	if(nrhs != 6) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. PI initial vector
	2. M transition matrix
	3. N mixture matrix 
	4. MU mean matrix
	5. SIGMA standard error matrix
	6. seq len

 ---   Note : Dimensions are already IN the model !!! ---
	*/






	// Here Just Determine Dimensions ! 
	in = mxGetPr(prhs[0]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the x dim
	in = mxGetPr(prhs[2]);   // Get the mixture matrix
	L = mxGetM(prhs[2]);
	y_dim = mxGetN(prhs[2]);     // Get the y dim


	hmm.miss_data = miss_data; // double .. 
	hmm.x_dim = x_dim;
	hmm.y_dim = y_dim;
	hmm.y_type = CONTINUOUS;
	hmm.cum_flag = 1;
	hmm.log_flag = 1;
	hmm.place_flag = place_flag; // currently no placing 
	hmm.use_bounds = use_bounds;	
	long special_models_flag = 0; 
	
	// Here init the model to get everything in the right place !!! 
	InitilizeModel(&hmm, x_dim, y_dim, CONTINUOUS, use_bounds, place_flag, special_models_flag, miss_data, 0);  // fold change flag is set to zero here,,
																			   // since we do not do learning anyway, 
																			   // and don't want to 'ruin' the mu's

	in = mxGetPr(prhs[0]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the x dim
  	for(i = 0; i < x_dim; i++)
		hmm.PI[i] = in[i];   // Copy the initial vector


	in = mxGetPr(prhs[1]);   // Get the transition matrix
	x_dim = mxGetM(prhs[1]);     // Get the x dim
	L = mxGetN(prhs[1]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");

  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix


	in = mxGetPr(prhs[2]);   // Get the mixture matrix
	L = mxGetM(prhs[2]);
	y_dim = mxGetN(prhs[2]);     // Get the y dim
	if(L != x_dim)
		mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[3]);   // Get the mean matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[4]);   // Get the std matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.SIGMA[i][j] = in[j*x_dim+i];   // Copy the std matrix


	in = mxGetPr(prhs[5]);   // Get the  sequence length
	seq_len = long(in[0]); 


	data.seq_len = seq_len; // this does in the data
	data.y_type = CONTINUOUS;
	data.miss_data = miss_data;
	data.dont_miss_prob = 0.5; // not important (only for simulations)

	ComputeModelAuxillaryParameters(&hmm); // Compute the cums ...

	// allocate memory
	for(i = 0; i < x_dim; i++)
	{
		data.x_vec = new long[seq_len];
		data.mix_vec = new long[seq_len];
		data.y_vec = new double[seq_len];

	}



	// Now Simulate The Sequence : 
	SimulateSequenceFromModel(&hmm, &data);

	// Now copy Data output vector  : 
	plhs[0] = mxCreateDoubleMatrix(data.seq_len, 1,  mxREAL);  // The Output vector 
	out = mxGetPr(plhs[0]); // output Viterbi
	for(i=0; i<seq_len; i++) 
		out[i] = data.y_vec[i];


	// New ! Now copy also the hidden vector !
	plhs[1] = mxCreateDoubleMatrix(data.seq_len, 1,  mxREAL);  // The Output vector 
	out = mxGetPr(plhs[1]); // output Viterbi
	for(i=0; i<seq_len; i++) 
		out[i] = data.x_vec[i];



  // free memory
	delete data.x_vec; 
	delete data.mix_vec; 
	delete data.y_vec; 





}
//#endif // MEX_TO_MATLAB







// Print everything about the model
long PrintModel2(hmm_model *hmm)
{
	long i;


	printf("HMM MODEL : \n");
	printf("X DIM : %ld\n", hmm->x_dim);


	if(hmm->y_type == DISCRETE)
	{
		printf("Y DIM : %ld\n", hmm->y_dim);
		printf("X is Discrete, Y is Discrete\n");
	}
	else
		printf("X is Discrete, Y is Mixture of %ld Gaussians\n", hmm->y_dim);


	printf("Initial X distribution : \n");
	PrintDoubleVec(hmm->PI, hmm->x_dim);

	printf("X Transision Matrix :\n");
	for(i = 0; i < hmm->x_dim; i++)
		PrintDoubleVec(hmm->M[i], hmm->x_dim);

	printf("X CUMCUM Transision Matrix :\n");
	for(i = 0; i < hmm->x_dim; i++)
		PrintDoubleVec(hmm->M_cum[i], hmm->x_dim);

	if(hmm->y_type == DISCRETE)
	{
		printf("X->Y Omission Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec(hmm->N[i], hmm->y_dim);


		printf("X->Y CUMCUM Omission Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec(hmm->N_cum[i], hmm->y_dim);
	}
	else
	{
		printf("X->Y Mixture Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec(hmm->N[i], hmm->y_dim);
		printf("Y Mean Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec(hmm->MU[i], hmm->y_dim);
		printf("Y Std. Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec(hmm->SIGMA[i], hmm->y_dim);
	}



	return 0;

}







