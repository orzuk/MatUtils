#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
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
/// Note : HMMs (hidden states) can be of different dimensions !!!!!!!  This doesn't change anything for the program
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
	int nRows, nCols, i, j, K, t, max_iters, num_starts;

	long seq_len, L; 
	long x_dim, y_dim, place_flag, miss_data, num_iters;
	long sum_x_dim, special_models_flag, p;
	double tolerance;
  
	double *in, *out;

	hmm_model hmm1;  // containts the first model
	hmm_model hmm2;  // containts the second model

	hmm_data data; // contains the data


	double KL_dist;

	 /* Check for proper number of arguments. Should be Fourteen 14 */
	if(nrhs != 14) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. Special Models Flag
	2. Length of sequence to simulate 
	3. Number of simulations to perform
	4. place_M matrices (only when special flag is on)
	5. First PI initial vector
	6. First M transition matrix
	7. First N mixture matrix 
	8. First MU mean matrix
	9. First SIGMA standard error matrix
	10. Second PI initial vector
	11. Second M transition matrix
	12. Second N mixture matrix 
	13. Second MU mean matrix
	14. Second SIGMA standard error matrix

 ---   Note : Dimensions are already IN the model !!! ---
	*/

///////////////////////////////////////////////////////////////
//// Read the data ...
///////////////////////////////////////////////////////////////
	in = mxGetPr(prhs[0]);		// Get the special models flag
	special_models_flag = long(in[0]);		// Get the special models flag	
	in = mxGetPr(prhs[1]);		// Get the expression vector length
	seq_len = long(in[0]);		// Get the expression vector length
	in = mxGetPr(prhs[2]);		// Get the number of iterations
	num_iters = long(in[0]);	// Get the number of iterations


	
	if(special_models_flag)
	{
		hmm1.x_dim2 = 2; hmm2.x_dim2 = 2;
		for(i = 0; i < hmm1.x_dim2; i++)		  // allocate memory
			for(j = 0; j < hmm1.x_dim2; j++)
			{
				hmm1.place_M[i][j] = new double[seq_len];
				hmm1.place_M_cum[i][j] = new double[seq_len]; 
				hmm2.place_M[i][j] = new double[seq_len];
				hmm2.place_M_cum[i][j] = new double[seq_len]; 
			}
		in = mxGetPr(prhs[3]);   // Get the place_M matrices vector
		for(p=0; p < seq_len; p++)  // copy the place_M matrices
			for(i = 0; i < hmm1.x_dim2; i++)
				for(j = 0; j < hmm1.x_dim2; j++)
				{
					hmm1.place_M[i][j][p] = in[p*hmm1.x_dim2*hmm1.x_dim2+j*hmm1.x_dim2+i];   // Copy the place_M matrix
					hmm2.place_M[i][j][p] = hmm1.place_M[i][j][p];  // Copy the place_M matrix
				}
	}



	in = mxGetPr(prhs[4]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[4]), mxGetN(prhs[4]));     // Get the x dim
	for(i = 0; i < x_dim; i++)
		hmm1.PI[i] = in[i];   // Copy the initial vector
	sum_x_dim = x_dim + (x_dim-1)*special_models_flag;

	in = mxGetPr(prhs[5]);   // Get the transition matrix
	x_dim = mxGetM(prhs[5]);     // Get the x dim
	L = mxGetN(prhs[5]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm1.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix

	in = mxGetPr(prhs[6]);   // Get the mixture matrix
	L = mxGetM(prhs[6]);
	y_dim = mxGetN(prhs[6]);     // Get the y dim
	if(L != sum_x_dim   )
		mexErrMsgTxt("Usage: eeeeeemission and transition matrices must have same number of rows ");
  	for(i = 0; i < sum_x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm1.N[i][j] = in[j*sum_x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[7]);   // Get the mean matrix
  	for(i = 0; i < sum_x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm1.MU[i][j] = in[j*sum_x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[8]);   // Get the std matrix
  	for(i = 0; i < sum_x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm1.SIGMA[i][j] = in[j*sum_x_dim+i];   // Copy the std matrix

	place_flag = special_models_flag; // space flag only in the special SNP case 
	miss_data = 0; 
	hmm1.seq_len = seq_len;
	hmm1.miss_data = miss_data; // double .. 
	hmm1.place_flag = place_flag; // copy the placing 
	hmm1.special_models_flag = special_models_flag;
	hmm1.x_dim = x_dim;
	hmm1.y_dim = y_dim;
	hmm1.y_type = CONTINUOUS;
	hmm1.cum_flag = 1;
	hmm1.log_flag = 1;


	

///////////////////////////////////////////////////////////////
//// Read the data for the second HMM ...
///////////////////////////////////////////////////////////////
	in = mxGetPr(prhs[9]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[9]), mxGetN(prhs[9]));     // Get the x dim
  	for(i = 0; i < x_dim; i++)
		hmm2.PI[i] = in[i];   // Copy the initial vector
	sum_x_dim = x_dim + (x_dim-1)*special_models_flag;

	in = mxGetPr(prhs[10]);   // Get the transition matrix
	x_dim = mxGetM(prhs[10]);     // Get the x dim
	L = mxGetN(prhs[10]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm2.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix

	in = mxGetPr(prhs[11]);   // Get the mixture matrix
	L = mxGetM(prhs[11]);
	y_dim = mxGetN(prhs[11]);     // Get the y dim
	if(L != sum_x_dim)
		mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
  	for(i = 0; i < sum_x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm2.N[i][j] = in[j*sum_x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[12]);   // Get the mean matrix
  	for(i = 0; i < sum_x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm2.MU[i][j] = in[j*sum_x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[13]);   // Get the std matrix
  	for(i = 0; i < sum_x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm2.SIGMA[i][j] = in[j*sum_x_dim+i];   // Copy the std matrix
///////////////////////////////////////////////////////////////


	data.seq_len = seq_len; // this does in the data
	data.y_type = CONTINUOUS;
	data.miss_data = miss_data;
	data.dont_miss_prob = 0.5; // not important (only for simulations)
	hmm2.seq_len = seq_len;
	hmm2.miss_data = miss_data; // double .. 
	hmm2.x_dim = x_dim;
	hmm2.y_dim = y_dim;
	hmm2.y_type = CONTINUOUS;
	hmm2.cum_flag = 1;
	hmm2.log_flag = 1;
	hmm2.place_flag = place_flag; // copy the placing 
	hmm2.special_models_flag = special_models_flag;


	hmm1.seq_len = seq_len;
	hmm1.miss_data = miss_data; // double .. 
	hmm1.place_flag = place_flag; // copy the placing 
	hmm1.special_models_flag = special_models_flag;
	hmm1.x_dim = x_dim;
	hmm1.y_dim = y_dim;
	hmm1.y_type = CONTINUOUS;
	hmm1.cum_flag = 1;
	hmm1.log_flag = 1;


	// Here init the model to get everything in the right place !!! 
	ComputeModelAuxillaryParameters(&hmm1); // Compute the cums ...
	ComputeModelAuxillaryParameters(&hmm2); // Compute the cums ...

	// Now Call the c function to do the job for you and give the distance ..


	if(hmm1.special_models_flag)
	{
		printf("Start SNPS KLD SP. Flags %ld %ld SeqLen %ld MISS DATA %ld XDIM %ld PI CUM\n", 
			hmm1.special_models_flag, hmm2.special_models_flag, data.seq_len, data.miss_data, hmm1.x_dim); fflush(stdout);
		for(i=0;i< hmm1.x_dim; i++)
			printf("%lf ", hmm1.PI_cum[i]);
		printf("\nSSSSS Tart\n");
		for(i=0;i< hmm1.x_dim; i++)
		{
			for(j=0;j< hmm1.x_dim; j++)
				printf("%lf ", hmm1.M_cum[i][j]);
			printf("\n");
		}

		KL_dist = 0.9999;
		
		ComputeKLDistanceSNPs(&hmm1, &hmm2, &data, num_iters, &KL_dist);
///	return; // early termination
	}
	else
	{

/***
		printf("Start Regular KLD SP. Flags %ld %ld SeqLen %ld MISS DATA %ld XDIM %ld PI CUM\n", 
			hmm1.special_models_flag, hmm2.special_models_flag, data.seq_len, data.miss_data, hmm1.x_dim); fflush(stdout);
		for(i=0;i< hmm1.x_dim; i++)
			printf("%lf ", hmm1.PI_cum[i]);

		printf("\nSSSSS Tart\n");
		for(i=0;i< hmm1.x_dim; i++)
		{
			for(j=0;j< hmm1.x_dim; j++)
				printf("%lf ", hmm1.M_cum[i][j]);
			printf("\n");
		}
**/
		ComputeKLDistance(&hmm1, &hmm2, &data, num_iters, &KL_dist);
	}
		printf("Ended  KLD\n"); fflush(stdout);


	// Now copy KL distance  output : 
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);  // The output number 
	out = mxGetPr(plhs[0]); // output number
	out[0] = KL_dist;

//////////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************ All Phases 



************************************************************************************ Ended All Phases ****/
//////////////////////////////////////////////////////////////////////////////////////////////////////


	if(special_models_flag)
		for(i = 0; i < hmm1.x_dim2; i++)		  // delete memory
			for(j = 0; j < hmm1.x_dim2; j++)
			{
				delete hmm1.place_M[i][j];
				delete hmm1.place_M_cum[i][j]; 
				delete hmm2.place_M[i][j];
				delete hmm2.place_M_cum[i][j]; 
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



