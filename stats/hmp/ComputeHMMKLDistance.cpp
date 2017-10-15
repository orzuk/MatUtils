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
long PrintDoubleVec2(double *vec, long len);

///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For HMM Chromosome project - 01.08.2004
/// Compute the distance between two HMMs 
/// Note : HMMs can be of different dimensions !!!!!!! 
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
	int nRows, nCols, i, j, K, t, max_iters, num_starts;

	long seq_len, L; 
	long x_dim, y_dim, place_flag, miss_data, num_iters, y_type;
	double tolerance;
  
	double *in, *out;

	hmm_model hmm1;  // containts the first model
	hmm_model hmm2;  // containts the second model


	hmm_data data; // contains the data


	double KL_dist;

	 /* Check for proper number of arguments. Should be Seven 7 */
	if(nrhs != 13) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
  	1. First PI initial vector
	2. First M transition matrix
	3. First N mixture matrix 
	4. First MU mean matrix
	5. First SIGMA standard error matrix
	6. Second PI initial vector
	7. Second M transition matrix
	8. Second N mixture matrix 
	9. Second MU mean matrix
	10. Second SIGMA standard error matrix
	11. Length of sequence to simulate 
	12. Number of simulations to perform
	13. New ! flag saying if we do Gaussians or Discrete ! 


 ---   Note : Dimensions are already IN the model !!! ---
	*/



//	printf("Start ...\n");
//	fflush(stdout); 


///////////////////////////////////////////////////////////////
//// Read the data ...
///////////////////////////////////////////////////////////////
	in = mxGetPr(prhs[0]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the x dim

//	printf("XDIM IS : %ld\n", x_dim); 
//	fflush(stdout); 

  	for(i = 0; i < x_dim; i++)
		hmm1.PI[i] = in[i];   // Copy the initial vector


	in = mxGetPr(prhs[1]);   // Get the transition matrix
	x_dim = mxGetM(prhs[1]);     // Get the x dim
	L = mxGetN(prhs[1]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");

  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm1.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix


	in = mxGetPr(prhs[2]);   // Get the mixture matrix
	L = mxGetM(prhs[2]);
	y_dim = mxGetN(prhs[2]);     // Get the y dim
	if(L != x_dim)
		mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm1.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[3]);   // Get the mean matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm1.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[4]);   // Get the std matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm1.SIGMA[i][j] = in[j*x_dim+i];   // Copy the std matrix


	// New : Get a flag saying if we are in the discrete or gaussian case. If discrete, mu and sigma have no meaning
	in = mxGetPr(prhs[12]);		// Get the discrete/gaussian flag
	y_type = long(in[0]);		// Get the discrete/gaussian flag


	place_flag = 0; 
	miss_data = 0; 

	hmm1.miss_data = miss_data; // double .. 
	hmm1.place_flag = place_flag; // currently no placing 
	hmm1.x_dim = x_dim;
	hmm1.y_dim = y_dim;
	hmm1.y_type = y_type; // CONTINUOUS;
	hmm1.cum_flag = 1;
	hmm1.log_flag = 1;


//	printf("Read 1st data ...\n");
///////////////////////////////////////////////////////////////
//// Read the data for the second HMM ...
///////////////////////////////////////////////////////////////
	in = mxGetPr(prhs[5]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[5]), mxGetN(prhs[5]));     // Get the x dim

  	for(i = 0; i < x_dim; i++)
		hmm2.PI[i] = in[i];   // Copy the initial vector


	in = mxGetPr(prhs[6]);   // Get the transition matrix
	x_dim = mxGetM(prhs[6]);     // Get the x dim
	L = mxGetN(prhs[6]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");

  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm2.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix


	in = mxGetPr(prhs[7]);   // Get the mixture matrix
	L = mxGetM(prhs[7]);
	y_dim = mxGetN(prhs[7]);     // Get the y dim
	if(L != x_dim)
		mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm2.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[8]);   // Get the mean matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm2.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[9]);   // Get the std matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm2.SIGMA[i][j] = in[j*x_dim+i];   // Copy the std matrix


///////////////////////////////////////////////////////////////

//	printf("Read 2nd data ...\n");


	in = mxGetPr(prhs[10]);		// Get the missing data flag
	seq_len = long(in[0]);		// Get the expression vector length

	in = mxGetPr(prhs[11]);		// Get the missing data flag
	num_iters = long(in[0]);		// Get the expression vector length

	


//	printf("seq len : %ld num iters : %ld \n", seq_len, num_iters);

//	data.y_vec = new double[seq_len];  // allocate the y vector 
//	data.y_vec_int = new long[seq_len];  // allocate the y integer vector 
//	data.x_vec = new long[seq_len];  // allocate the x vector
//	data.loc_vec = new long[seq_len];  // allocate locations vector
//	data.loc_diff_vec = new long[seq_len];	// allocate diff location vectors
//	data.mix_vec = new long[seq_len];  // allocate the mix vec 

	


	data.seq_len = seq_len; // this does in the data
	data.y_type = y_type; // CONTINUOUS;
	data.miss_data = miss_data;
	data.dont_miss_prob = 0.5; // not important (only for simulations)



	hmm2.miss_data = miss_data; // double .. 
	hmm2.x_dim = x_dim;
	hmm2.y_dim = y_dim;
	hmm2.y_type = y_type; // CONTINUOUS;
	hmm2.cum_flag = 1;
	hmm2.log_flag = 1;
	hmm2.place_flag = place_flag; // currently no placing 




	// Here init the model to get everything in the right place !!! 
	


	printf("START AUX\n");

	ComputeModelAuxillaryParameters(&hmm1); // Compute the cums ...
	ComputeModelAuxillaryParameters(&hmm2); // Compute the cums ...


//	printf("HMM1 : \n");
//	PrintModel2(&hmm1);
//	printf("HMM2 : \n");
//	PrintModel2(&hmm2);

//	printf("Calling KLKKLKLK ...\n");
	
	// Now Call the c function to do the job for you and give the distance ..
	ComputeKLDistance(&hmm1, &hmm2, &data, num_iters, &KL_dist); // Where is this function?? 

	printf("Finished  KLKKLKLK Got Distance : %lf ...\n", KL_dist);


	// Now copy KL distance  output : 
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);  // The output number 
	out = mxGetPr(plhs[0]); // output number
	out[0] = KL_dist;

	printf("Now deleting ...\n");


/***
	plhs[0] = mxCreateDoubleMatrix(num_seqs, seqs_len-L+1, mxREAL);
	out = mxGetPr(plhs[0]); // output the scores array 
	for(i=0; i<num_seqs; i++) 
		for (j=0; j<seqs_len-L+1; j++)
			out[j*num_seqs+i] = scores[i][j]; // copy scores array to output

***/		

  // free memory

//	delete data.x_vec;
//	delete data.y_vec;
//	delete data.y_vec_int;
//	delete data.loc_vec;
//	delete data.loc_diff_vec;
//	delete data.mix_vec;
	
		

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
	PrintDoubleVec2(hmm->PI, hmm->x_dim);

	printf("X Transision Matrix :\n");
	for(i = 0; i < hmm->x_dim; i++)
		PrintDoubleVec2(hmm->M[i], hmm->x_dim);

	printf("X CUMCUM Transision Matrix :\n");
	for(i = 0; i < hmm->x_dim; i++)
		PrintDoubleVec2(hmm->M_cum[i], hmm->x_dim);

	if(hmm->y_type == DISCRETE)
	{
		printf("X->Y Omission Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec2(hmm->N[i], hmm->y_dim);


		printf("X->Y CUMCUM Omission Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec2(hmm->N_cum[i], hmm->y_dim);
	}
	else
	{
		printf("X->Y Mixture Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec2(hmm->N[i], hmm->y_dim);
		printf("Y Mean Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec2(hmm->MU[i], hmm->y_dim);
		printf("Y Std. Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec2(hmm->SIGMA[i], hmm->y_dim);
	}



	return 0;

}

long PrintDoubleVec2(double *vec, long len)
{
	long i;
	for(i = 0; i < len; i++)
	{
		printf("%lf ", vec[i]);
		if((i%PRINT_IN_ROW) == PRINT_IN_ROW-1)
			printf("\n");
	}

	printf("\n");

	return 0; 
}




