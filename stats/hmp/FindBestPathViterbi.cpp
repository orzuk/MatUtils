#include <stdlib.h>
#include <stdio.h>
#include <iostream> // new convension: no .h
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
	long x_dim, y_dim, place_flag, miss_data;
	double tolerance;
  
	double *in, *out;

//	long *y_vec_int; // dummy !!! 
//	double *y_vec; 
//	double *loc_vec;  // location on chromosome vector - Currently not used !!!! 
	
	hmm_model hmm;  // containts the model

	hmm_data data; // contains the data

	double *y_cond_tabs[MAX_X_VALS]; // helpful conditional tables
	double *alpha[MAX_X_VALS];
	double *beta[MAX_X_VALS];
	double *gamma[MAX_X_VALS];
	double *gammagamma[MAX_X_VALS][MAX_Y_VALS];
	double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS];


	long *opt_x_vec_out;  // optimal viterbi path

	double *scale, *scale_cum, *scale_exp;
	double cur_score;

	 /* Check for proper number of arguments. Should be Eleven 11 */
	if(nrhs != 11) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. Expression vector
	2. Chromosomal location vector
	3. missing data flag 
	4. PI initial vector
	5. M transition matrix
	6. N mixture matrix 
	7. MU mean matrix
	8. SIGMA standard error matrix
	9. place flag
	10. place means
	11. place sigmas

 ---   Note : Dimensions are already IN the model !!! ---
	*/



	in = mxGetPr(prhs[0]);   // Get the expression vector
	seq_len = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the expression vector length

	data.y_vec = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy expression vector
		data.y_vec[t] = in[t];  
//	data.y_vec_int = new long[seq_len];
//	for(t = 0; t < seq_len; t++)   // copy expression to integer vector - maybe this will help (?)
//		data.y_vec_int[t] = long(in[t]);  





	in = mxGetPr(prhs[1]);   // Get the location vector
	L  = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));    // Get the number of sequences (genes)
	if(L != seq_len)
		mexErrMsgTxt("Usage: Exp. and Loc. vector length must be the same");

	
	data.loc_vec = new long[seq_len];
	for(t = 0; t < seq_len; t++)   // copy location vector
		data.loc_vec[t] = long(in[t]/NUM_BASEPAIRS_PER_HMM_TRANS);  

	data.loc_diff_vec = new long[seq_len];
	data.loc_diff_vec[0] = 0; 
	for(t = 1; t < seq_len; t++)   // compute location diff vector
		data.loc_diff_vec[t] = MIN(MAX(data.loc_vec[t]-data.loc_vec[t-1]-1, 0), MAX_POWER_OF_TRANS_MATRIX-2);  


	in = mxGetPr(prhs[2]);   // Get the missing data flag
	miss_data = long(in[0]);  


	in = mxGetPr(prhs[3]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[3]), mxGetN(prhs[3]));     // Get the x dim

  	for(i = 0; i < x_dim; i++)
		hmm.PI[i] = in[i];   // Copy the initial vector


	in = mxGetPr(prhs[4]);   // Get the transition matrix
	x_dim = mxGetM(prhs[4]);     // Get the x dim
	L = mxGetN(prhs[4]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");

  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix


	in = mxGetPr(prhs[5]);   // Get the mixture matrix
	L = mxGetM(prhs[5]);
	y_dim = mxGetN(prhs[5]);     // Get the y dim
	if(L != x_dim)
		mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[6]);   // Get the mean matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[7]);   // Get the std matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.SIGMA[i][j] = in[j*x_dim+i];   // Copy the std matrix

///	printf("xdim %ld ydim %ld\n", x_dim, y_dim);
///	printf("START WITH MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);


	in = mxGetPr(prhs[8]);   // Get the place flag
	place_flag = long(in[0]);


	if(place_flag)
	{
		// allocate memory 
		for(i = 0; i < x_dim; i++)		
			for(j = 0; j < y_dim; j++)
			{
				hmm.place_gaussian_mu[i][j] = new double[seq_len];
				hmm.place_gaussian_sigma[i][j] = new double[seq_len];
			}



		// Now copy
		in = mxGetPr(prhs[9]);   // Get the place means
		long arr_size = mxGetM(prhs[9]) * mxGetN(prhs[9]);
		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
			{
				for(t = 0; t < seq_len; t++)
				{
			//		if(t+seq_len*j+seq_len*y_dim*i >= arr_size)
			//			printf("VITTTTERROR !! TOO LARGE INDEX %ld\n", t+seq_len*j+seq_len*y_dim*i);
					hmm.place_gaussian_mu[i][j][t] = in[t+seq_len*j+seq_len*y_dim*i]; // Copy the place means
				}
			}
		// Now copy
		in = mxGetPr(prhs[10]);   // Get the place sigmas
		for(i = 0; i < x_dim; i++)	
			for(j = 0; j < y_dim; j++)
				for(t = 0; t < seq_len; t++)
					hmm.place_gaussian_sigma[i][j][t] = in[t+seq_len*j+seq_len*y_dim*i]; // Copy the place sigmas


	}


	data.seq_len = seq_len; // this does in the data
	data.y_type = CONTINUOUS;
	data.miss_data = miss_data;
	data.dont_miss_prob = 0.5; // not important (only for simulations)


	hmm.miss_data = miss_data; // double .. 
	hmm.x_dim = x_dim;
	hmm.y_dim = y_dim;
	hmm.y_type = CONTINUOUS;
	hmm.cum_flag = 1;
	hmm.log_flag = 1;
	hmm.place_flag = place_flag; // currently no placing 
	long special_models_flag = 0; 

///	printf("BEFORE INIT MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);



	// Here init the model to get everything in the right place !!! (use_bound is set to 0 since we dont need them in Viterbi)
	InitilizeModel(&hmm, x_dim, y_dim, CONTINUOUS, 0, place_flag, special_models_flag, miss_data, 0);  // fold change flag is set to zero here,,
																			   // since we do not do learning anyway, 
																			   // and don't want to 'ruin' the mus
	

	

	// Now write again the parameters !!!!!
	in = mxGetPr(prhs[3]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[3]), mxGetN(prhs[3]));     // Get the x dim

  	for(i = 0; i < x_dim; i++)
		hmm.PI[i] = in[i];   // Copy the initial vector


	in = mxGetPr(prhs[4]);   // Get the transition matrix
	x_dim = mxGetM(prhs[4]);     // Get the x dim
	L = mxGetN(prhs[4]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");

  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix


	in = mxGetPr(prhs[5]);   // Get the mixture matrix
	L = mxGetM(prhs[5]);
	y_dim = mxGetN(prhs[5]);     // Get the y dim
	if(L != x_dim)
		mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[6]);   // Get the mean matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[7]);   // Get the std matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.SIGMA[i][j] = in[j*x_dim+i];   // Copy the std matrix




///	printf("AFTER INIT MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);	
	ComputeModelAuxillaryParameters(&hmm); // Compute the cums ...


////	printf("AFTER AUXIL MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);
	// allocate memory
	for(i = 0; i < x_dim; i++)
	{
		y_cond_tabs[i] = new double[seq_len];
		alpha[i] = new double[seq_len];
		beta[i] = new double[seq_len];
		gamma[i] = new double[seq_len];

		for(j = 0; j < y_dim; j++)
		{
			gammagamma[i][j] = new double[seq_len];
			y_mu_square_exp_vecs[i][j] = new double[seq_len];
		}
	}
	opt_x_vec_out = new long[seq_len];
	scale = new double[seq_len];
	scale_cum = new double[seq_len];
	scale_exp = new double[seq_len];

	// Call many c functions to do the job for you ..

	Compute_y_cond_tabs(&hmm, &data, y_cond_tabs, y_mu_square_exp_vecs);  // OK


	// Now find the probs
	forward(&hmm, &data, y_cond_tabs, alpha, scale, scale_cum, scale_exp, &cur_score);  // PROBLEM FIRST ALPHA !!! 
	backward(&hmm, &data, scale_exp, y_cond_tabs, beta);		

	ComputeMarginalGamma(&hmm, &data, y_mu_square_exp_vecs, alpha, beta, scale_exp, y_cond_tabs,
						 gamma, gammagamma);


	// Print the Viterbi output vector
///	printf("Marginal Out : \n");
///	for(i=0; i<seq_len; i++)
///		printf("%lf ,", gamma[1][i]);
///	printf("\n");
///	printf("NOW MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);



	double *x_loglike, *x_place_loglike, *y_loglike; 
	Viterbi( &hmm, &data, y_cond_tabs, y_cond_tabs, opt_x_vec_out, x_loglike, x_place_loglike, y_loglike);

	// Print the Viterbi output vector
/***
	printf("Vit Out : \n");
	for(i=0; i<seq_len; i++)
		printf("%ld ,", opt_x_vec_out[i]);
	printf("\n");
***/


	// Now copy Viterbi output vector  : 
	plhs[0] = mxCreateDoubleMatrix(1, seq_len, mxREAL);  // The Viterbi vector 
	out = mxGetPr(plhs[0]); // output Viterbi
	for(i=0; i<seq_len; i++) 
		out[i] = opt_x_vec_out[i];


	// Now copy gamma probs. output vector
	plhs[1] = mxCreateDoubleMatrix(x_dim, seq_len, mxREAL);  // The gamma vector 
	out = mxGetPr(plhs[1]); // output gamma
	for(i = 0;i < x_dim; i++)
		for(t=0; t<seq_len; t++) 
			out[i+t*x_dim] = gamma[i][t];






/***
	plhs[0] = mxCreateDoubleMatrix(num_seqs, seqs_len-L+1, mxREAL);
	out = mxGetPr(plhs[0]); // output the scores array 
	for(i=0; i<num_seqs; i++) 
		for (j=0; j<seqs_len-L+1; j++)
			out[j*num_seqs+i] = scores[i][j]; // copy scores array to output

***/		

  // free memory

	for(i = 0; i < x_dim; i++)
	{
		delete y_cond_tabs[i];
		delete alpha[i];
		delete beta[i];
		delete gamma[i];
		for(j = 0; j < y_dim; j++)
		{
			delete gammagamma[i][j];
			delete y_mu_square_exp_vecs[i][j];
		}
	}


	delete data.y_vec;
//	delete data.y_vec_int;
	delete data.loc_vec;
	delete data.loc_diff_vec;

	
	delete scale;
	delete opt_x_vec_out;
	delete scale_cum; 
	delete scale_exp;
		

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







