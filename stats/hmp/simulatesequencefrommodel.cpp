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

///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For HMM Chromosome project - 11.07.2004
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
	int nRows, nCols, i, j, p, K, t, max_iters, num_starts;

	long seq_len, L; 
	long x_dim, y_dim; // , place_flag, miss_data;
	double tolerance;
  
	double *in, *out;

//	long *y_vec_int; // dummy !!! 
//	double *y_vec; 
//	double *loc_vec;  // location on chromosome vector - Currently not used !!!! 
	
	hmm_model hmm;  // containts the model

	hmm_data data; // contains the data


	long *opt_x_vec_out;  // optimal viterbi path

	double *scale, *scale_cum, *scale_exp;
	double cur_score;
	
	long miss_data = 0; // Don't miss data for now !!!! 
	long place_flag, special_models_flag; // Use different transition matrices ?
	long use_bounds = 0; // don't use bounds for now !!! 


	
	/* Check for proper number of arguments. Should be Nine 9 */
	if(nrhs != 9) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. PI initial vector
	2. M transition matrix
	3. N mixture matrix 
	4. MU mean matrix
	5. SIGMA standard error matrix
	6. place flag
	7. place M transition matrices
	8. special models flag
	9. seq len

 ---   Note : Dimensions are already IN the model !!! ---
	*/


	// Here Just Determine Dimensions ! 
	in = mxGetPr(prhs[0]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the x dim
	in = mxGetPr(prhs[2]);   // Get the mixture matrix
	L = mxGetM(prhs[2]);
	y_dim = mxGetN(prhs[2]);     // Get the y dim

	in = mxGetPr(prhs[5]);   // Get the  place flag
	place_flag = long(in[0]); 

	in = mxGetPr(prhs[7]);   // Get the  special flag
	special_models_flag = long(in[0]); 

	hmm.miss_data = miss_data; // double .. 
	hmm.x_dim = x_dim;
	hmm.y_dim = y_dim;
	hmm.y_type = CONTINUOUS;
	hmm.cum_flag = 1;
	hmm.log_flag = 1;
	hmm.place_flag = place_flag;  
	hmm.use_bounds = use_bounds;
	hmm.special_models_flag = special_models_flag;

	
	in = mxGetPr(prhs[8]);   // Get the  sequence length
	seq_len = long(in[0]); 
	hmm.seq_len = seq_len;
//	printf("Start ALLOCATING\n");
	// allocate memory for mu and sigma 
	if(hmm.special_models_flag == SNPS_ONE_SAMPLE)
	{
		hmm.x_dim2 = 2; // set the binary genotype variables 
		for(i = 0; i < hmm.x_dim2; i++)
			for(j = 0; j < hmm.x_dim2; j++)
			{
				hmm.place_M[i][j] = new double[seq_len];
				hmm.place_M_cum[i][j] = new double[seq_len];
			}
	}
//	printf("Start Copying\n");

	in = mxGetPr(prhs[6]);  // Get the place transition matrices
	if(hmm.special_models_flag == SNPS_ONE_SAMPLE)
		for(p=0; p < seq_len; p++)
			for(i = 0; i < hmm.x_dim2; i++)
				for(j = 0; j < hmm.x_dim2; j++)
					hmm.place_M[i][j][p] = in[p*hmm.x_dim2*hmm.x_dim2+j*hmm.x_dim2+i];   // Copy the place_M matrix


//	printf("THE M MATRIX : %lf %lf %lf %lf\n", hmm.place_M[0][0][0], hmm.place_M[0][1][0], hmm.place_M[1][0][0], hmm.place_M[1][1][0]); 

	// Here init the model to get everything in the right place !!! 				
	InitilizeModel(&hmm, x_dim, y_dim, CONTINUOUS, use_bounds, place_flag, special_models_flag, miss_data, 0);  // fold change flag is set to zero here,,																			   // since we do not do learning anyway, 
																			   // and don't want to 'ruin' the mus

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

	
	if(hmm.special_models_flag == SNPS_ONE_SAMPLE)
	{

		if(L != 2*x_dim-1)
		{
			printf("dim: %d 2*dim-1: %d L %d\n", x_dim, 2*x_dim-1, L);
			mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows 2*dim-1");
		}
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.N[i][j] = in[j*(2*x_dim-1)+i];   // Copy the mixture matrix

		in = mxGetPr(prhs[3]);   // Get the mean matrix
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*(2*x_dim-1)+i];   // Copy the mean matrix
	
		in = mxGetPr(prhs[4]);   // Get the std matrix
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = in[j*(2*x_dim-1)+i];   // Copy the std matrix
	}
	else
	{
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
	}
	data.seq_len = seq_len; // this does in the data
	data.y_type = CONTINUOUS;
	data.miss_data = miss_data;
	data.dont_miss_prob = 0.5; // not important (only for simulations)

//	printf("First the Moments: mu %lf %lf %lf sigma %lf %lf %lf\n", hmm.MU[1][0], hmm.MU[1][0], hmm.MU[2][0], 
//														  hmm.SIGMA[0][0], hmm.SIGMA[1][0], hmm.SIGMA[2][0]);


	ComputeModelAuxillaryParameters(&hmm); // Compute the cums ...

	// allocate memory
	data.x_vec = new long[seq_len];
	data.mix_vec = new long[seq_len];
	data.y_vec = new double[seq_len];
	data.y_vecB = new double[seq_len]; // This is for the special SNP case


/**	printf("Start simulating\n");
	printf("Dimensions: x_dim %d x_dim2 %d y_dim %d\n", hmm.x_dim, hmm.x_dim2, hmm.y_dim);
	printf("Moments: mu %lf %lf %lf sigma %lf %lf %lf\n", hmm.MU[0][0], hmm.MU[1][0], hmm.MU[2][0], 
														  hmm.SIGMA[0][0], hmm.SIGMA[1][0], hmm.SIGMA[2][0]);
/***

	printf("The place cumsum\n");
	for(i=0;i<100;i++)
		printf("%lf %lf %lf %lf | %lf %lf %lf %lf | SUM %lf \n", 
		hmm.place_M[0][0][i], hmm.place_M[0][1][i], 
		hmm.place_M[1][0][i], hmm.place_M[1][1][i], 
		hmm.place_M_cum[0][0][i], hmm.place_M_cum[0][1][i], 
		hmm.place_M_cum[1][0][i], hmm.place_M_cum[1][1][i], 
		hmm.place_M[0][0][i] + hmm.place_M[0][1][i] + 
		hmm.place_M[1][0][i] + hmm.place_M[1][1][i]);

	printf("\nSpecial Models Flag : %ld  Miss data flag %ld Place Flag %ld DIM2 IS %ld\n", 
		hmm.special_models_flag, data.miss_data, hmm.place_flag, hmm.x_dim2); 
/***/

/***	
	// Now Simulate The Sequence : 
	printf("Before SIM, PI csumsum:\n");
	for(i=0; i < hmm.x_dim; i++)
		printf("%lf ", hmm.PI_cum[i]);
	printf("\n");
	printf("Before SIM, N csumsum:\n");
	for(i=0; i < (2*hmm.x_dim-1); i++)
		printf("%lf ", hmm.N_cum[i][0]);
	printf("\n");
***/


	SimulateSequenceFromModel(&hmm, &data);

	/***
	printf("THE FLAGS: SPECIAL %ld PLACE %ld\n", hmm.special_models_flag, hmm.place_flag);
	printf("InsideSIM Data out X-VEC\n");
	for(i=0; i<100; i++) 
		printf("%ld, ", data.x_vec[i]);
	printf("\nData out Y_VEC:\n");
	for(i=0; i<100; i++) 
		printf("%lf, ", data.y_vec[i]);
	printf("\nData out Y_VEC B:\n");
	for(i=0; i<100; i++) 
		printf("%lf, ", data.y_vecB[i]);
	printf("\n\n");
/***/

	// New ! Now copy also the hidden vector !
	plhs[0] = mxCreateDoubleMatrix(data.seq_len, 1,  mxREAL);  // The X Output vector 
	out = mxGetPr(plhs[0]); // output Viterbi
	for(i=0; i<seq_len; i++) 
		out[i] = data.x_vec[i];

	// Now copy Data output vector  : 
	plhs[1] = mxCreateDoubleMatrix(data.seq_len, 1,  mxREAL);  // The Y_A Output vector 
	out = mxGetPr(plhs[1]); // output Observations
	for(i=0; i<seq_len; i++) 
		out[i] = data.y_vec[i];

	// Now copy the second Data output vector  : 
	plhs[2] = mxCreateDoubleMatrix(data.seq_len, 1,  mxREAL);  // The Y_B Output vector 
	out = mxGetPr(plhs[2]); // output Observations
	if(hmm.special_models_flag == SNPS_ONE_SAMPLE)
		for(i=0; i<seq_len; i++) 
			out[i] = data.y_vecB[i];
	else // copy the same data twice
		for(i=0; i<seq_len; i++) 
			out[i] = data.y_vec[i];


  // free memory
	delete data.x_vec; 
	delete data.mix_vec; 
	delete data.y_vec; 
	delete data.y_vecB; 
	if(hmm.special_models_flag == SNPS_ONE_SAMPLE)
	{
		for(i = 0; i < hmm.x_dim2; i++)
			for(j = 0; j < hmm.x_dim2; j++)
			{
				delete hmm.place_M[i][j];
				delete hmm.place_M_cum[i][j];
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







