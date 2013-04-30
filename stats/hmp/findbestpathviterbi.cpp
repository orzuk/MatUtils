#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "general.h"
#include "hmm_chrom_funcs.h"


// #define DEBUG
#undef DEBUG

long PrintModel2(hmm_model *hmm);

///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For HMM Chromosome project - 11.07.2004
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  

	// Open again and start reading
#ifdef DEBUG
	double start_time, end_time;
	start_time = clock();
	char *debug_p = "DEBUG.txt";
	FILE *debug_f;
	debug_f = fopen(debug_p, "w");
	fprintf(debug_f, "Starting\n");
	fclose(debug_f);
#endif

	long nRows, nCols, i, j, K, t, max_iters, num_starts;
	long i1,i2,j1,j2,i_geno,j_geno,total_i_index;

	long seq_len, L; 
	long x_dim, x_dim2, y_dim, place_flag, miss_data, special_models_flag, snp_specific_flag;
	double tolerance;

	long num_gaussians, num_mix_gaussians_mu, num_mix_gaussians_sigma;

  
	double *in, *out;

//	double *loc_vec;  // location on chromosome vector - Currently not used !!!! 
	
	hmm_model hmm;  // containts the model
	hmm_data data; // contains the data

	double *y_cond_tabs[MAX_2D_X_VALS]; // helpful conditional tables
	double *y_cond_tabsB[MAX_2D_X_VALS]; // helpful conditional tables
	double *alpha[36];
	double *beta[36];
	double *gamma[36];
	double *phi[36][MAX_Y_VALS];
	double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS];
	double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS];



	long *opt_x_vec_out;  // optimal viterbi path

	double *scale, *scale_cum, *scale_exp;
	double cur_score;

	 /* Check for proper number of arguments. Should be Fourteen 14 */
	if(nrhs != 16) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. Expression vector
	2. Second Expression vector
	3. Chromosomal location vector
	4. missing data flag 
	5. PI initial vector
	6. M transition matrix
	7. N mixture matrix 
	8. MU mean matrix
	9. SIGMA standard error matrix
	10. place flag
	11. special model flag (NEW)	
	12. snp specific flag (NEW)
	13. place means
	14. place sigmas
	15. place sigmas invs
	16. place M transition matrices

 ---   Note : Dimensions are already IN the model !!! ---
	*/

//	printf("Starting Viteri\n"); fflush(stdout);

	in = mxGetPr(prhs[0]);   // Get the expression vector
	seq_len = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the expression vector length
	data.y_vec = new double[seq_len];
	if(data.y_vec == NULL) 
	{
		printf("Error: Failed to allocate data.y_vec\n");
		exit(0);
	}

	for(t = 0; t < seq_len; t++)   // copy expression vector
		data.y_vec[t] = in[t];  

	in = mxGetPr(prhs[1]);   // Get the additional expression vector if needed
	data.y_vecB = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy expression vector
		data.y_vecB[t] = in[t];  
		
	//	data.y_vec_int = new long[seq_len];
//	for(t = 0; t < seq_len; t++)   // copy expression to integer vector - maybe this will help (?)
//		data.y_vec_int[t] = long(in[t]);  



	in = mxGetPr(prhs[2]);   // Get the location vector
	L  = MAX(mxGetM(prhs[2]), mxGetN(prhs[2]));    // Get the number of sequences (genes)
	if(L != seq_len)
		mexErrMsgTxt("Usage: Exp. and Loc. vector length must be the same");

	
	data.loc_vec = new long[seq_len];
	for(t = 0; t < seq_len; t++)   // copy location vector
		data.loc_vec[t] = long(in[t]/NUM_BASEPAIRS_PER_HMM_TRANS);  

	data.loc_diff_vec = new long[seq_len];
	data.loc_diff_vec[0] = 0; 
	for(t = 1; t < seq_len; t++)   // compute location diff vector
		data.loc_diff_vec[t] = MIN(MAX(data.loc_vec[t]-data.loc_vec[t-1]-1, 0), MAX_POWER_OF_TRANS_MATRIX-2);  


	in = mxGetPr(prhs[3]);   // Get the missing data flag
	miss_data = long(in[0]);  


	in = mxGetPr(prhs[4]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[4]), mxGetN(prhs[4]));     // Get the x dim

  	for(i = 0; i < x_dim; i++)
		hmm.PI[i] = in[i];   // Copy the initial vector


	in = mxGetPr(prhs[5]);   // Get the transition matrix
	x_dim = mxGetM(prhs[5]);     // Get the x dim
	L = mxGetN(prhs[5]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix


	in = mxGetPr(prhs[9]);   // Get the place flag
	place_flag = long(in[0]);
	in = mxGetPr(prhs[10]);   // Get the special models flag
	special_models_flag = long(in[0]);
	in = mxGetPr(prhs[11]);   // Get the snp specific flag
	snp_specific_flag = long(in[0]);

	in = mxGetPr(prhs[6]);   // Get the mixture matrix
	L = mxGetM(prhs[6]);
	y_dim = mxGetN(prhs[6]);     // Get the y dim
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
		if(L != 2*x_dim-1)
		{
			printf("dim: %d 2*dim-1: %d L %d\n", x_dim, 2*x_dim-1, L);
			mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
		}
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.N[i][j] = in[j*(2*x_dim-1)+i];   // Copy the mixture matrix

		in = mxGetPr(prhs[7]);   // Get the mean matrix
	  	for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*(2*x_dim-1)+i];   // Copy the mean matrix

		in = mxGetPr(prhs[8]);   // Get the std matrix
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = MAX(MIN_SIGMA, in[j*(2*x_dim-1)+i]);   // Copy the st.d. matrix

//			for(i = 0; i < 2*x_dim-1; i++)
//			printf("MU %lf SIGMA %lf\n", hmm.MU[i][0], hmm.SIGMA[i][0]);
	
	}
	else
	{
		if(L != x_dim)
		{
			printf("WTFFFFF\n");
			mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
		}
  		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix
		in = mxGetPr(prhs[7]);   // Get the mean matrix
	  	for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix
		in = mxGetPr(prhs[8]);   // Get the std matrix
  		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = MAX(MIN_SIGMA, in[j*x_dim+i]);   // Copy the std matrix. Avoid too little std.s.
	}


///	printf("START WITH MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached 10\n");
	fclose(debug_f);
#endif

	if(place_flag)
	{	
		if(special_models_flag == SNPS_ONE_SAMPLE)
		{
			hmm.x_dim2 = 2; x_dim2=2; // set the binary genotype variables 		
			// allocate memory 
			for(i = 0; i < x_dim2; i++)		
				for(j = 0; j < x_dim2; j++)
				{
					hmm.place_M[i][j] = new double[seq_len];
					hmm.place_M_cum[i][j] = new double[seq_len];
				}

			// Now copy
			in = mxGetPr(prhs[15]);   // Get the place M transition matrices
			////long arr_size = mxGetM(prhs[12]) * mxGetN(prhs[12]);
			for(i = 0; i < x_dim2; i++)
				for(j = 0; j < x_dim2; j++)
					for(t = 0; t < seq_len-1; t++)
						hmm.place_M[i][j][t] = in[t*x_dim2*x_dim2+j*x_dim2+i]; // in[i+x_dim2*j+x_dim2*x_dim2*t]; // Copy the place M

 /***/


/***
			printf("THE PLACE_M MATRIX : %lf %lf %lf %lf\n", hmm.place_M[0][0][0], hmm.place_M[0][1][0], 
													   hmm.place_M[1][0][0], hmm.place_M[1][1][0]); 
			printf("THE PLACE_M MATRIX 10: %lf %lf %lf %lf\n", hmm.place_M[0][0][10], hmm.place_M[0][1][10], 
													   hmm.place_M[1][0][10], hmm.place_M[1][1][10]); 
			printf("THE PLACE_M MATRIX Last: %lf %lf %lf %lf\n", hmm.place_M[0][0][seq_len-2], hmm.place_M[0][1][seq_len-2], 
													   hmm.place_M[1][0][seq_len-2], hmm.place_M[1][1][seq_len-2]); 

			
			for(j=0; j < 20; j++)
				printf("THE PLACE_M %ld MATRIX : %lf %lf %lf %lf\n", j, hmm.place_M[0][0][j], hmm.place_M[0][1][j], 
													   hmm.place_M[1][0][j], hmm.place_M[1][1][j]); 

			
			/***/
		}
		
		if(snp_specific_flag || (!special_models_flag))
		{
			in = mxGetPr(prhs[12]);   // Get the place means
			if(snp_specific_flag)
			{
				hmm.gauss_dim = 2;
				num_gaussians = mxGetN(prhs[12])/2;    // Get the number of gaussians 
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "We've got %ld Gaussians\n", num_gaussians);
	fclose(debug_f);
#endif
			
				num_mix_gaussians_mu = 2;
				num_mix_gaussians_sigma = 3;
			}
			else
			{
				hmm.gauss_dim = 1;
				num_gaussians = x_dim;
				num_mix_gaussians_mu = y_dim;
				num_mix_gaussians_sigma = y_dim;
			}


			// allocate memory 
			for(i = 0; i < num_gaussians; i++)		
			{
				for(j = 0; j < num_mix_gaussians_mu; j++)
				{
					hmm.place_gaussian_mu[i][j] = new double[seq_len];
					if(hmm.place_gaussian_mu[i][j] == NULL) 
					{
						printf("Error: failed to allocate hmm.place_gaussian_mu\n");
						exit(99);
					}
				}
				for(j = 0; j < num_mix_gaussians_sigma; j++)
				{
					hmm.place_gaussian_sigma[i][j] = new double[seq_len];
					hmm.place_gaussian_sigma_inv[i][j] = new double[seq_len];

					if( (hmm.place_gaussian_sigma[i][j] == NULL) || (hmm.place_gaussian_sigma_inv[i][j] == NULL) )
					{
						printf("Error: failed to allocate hmm.place_gaussian_sigma\n");
						exit(99);
					}
				}
			}

			// Now copy
			in = mxGetPr(prhs[12]);   // Get the place means
			for(i = 0; i < num_gaussians; i++)
				for(j = 0; j < num_mix_gaussians_mu; j++)
					for(t = 0; t < seq_len; t++)
						hmm.place_gaussian_mu[i][j][t] = in[t+seq_len*i+seq_len*num_gaussians*j]; // Copy the place means
			// Now copy
			in = mxGetPr(prhs[13]);   // Get the place sigmas
			for(i = 0; i < num_gaussians; i++)	
				for(j = 0; j < num_mix_gaussians_sigma; j++)
					for(t = 0; t < seq_len; t++)
						hmm.place_gaussian_sigma[i][j][t] = in[t+seq_len*i+seq_len*num_gaussians*j]; // Copy the place sigmas
			// Now copy
			in = mxGetPr(prhs[14]);   // Get the place sigmas invs
			for(i = 0; i < num_gaussians; i++)	
				for(j = 0; j < num_mix_gaussians_sigma; j++)
					for(t = 0; t < seq_len; t++)
						hmm.place_gaussian_sigma_inv[i][j][t] = in[t+seq_len*i+seq_len*num_gaussians*j]; // Copy the place sigmas inverse
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	t=seq_len-1; j=num_mix_gaussians_sigma-1; i=num_gaussians-1;
	fprintf(debug_f, "last in index : %ld\n", t+seq_len*i+seq_len*num_gaussians*j);
	long dim_dim = mxGetM(prhs[14]) * mxGetN(prhs[14]);     // Get the array size  
	fprintf(debug_f, "In array size %ld\n", dim_dim);
	fclose(debug_f);
#endif



#ifdef DEBUG		
			// Print some of the first mu's and sigma's 			
			printf("\nGaussians MU SNP %ld:\n", t);
			for(j = 0; j < num_mix_gaussians_mu; j++)
				for(i = 0; i < num_gaussians; i++)
				{
					for(t = 0; t < seq_len/*10*/; t++)
						printf("%lf ", hmm.place_gaussian_mu[i][j][t]);
					printf("\n");
				}

			printf("\nGaussians SIGMA SNP %ld:\n", t);
			for(j = 0; j < num_mix_gaussians_sigma; j++)
				for(i = 0; i < num_gaussians; i++)
				{
					for(t = 0; t < seq_len/*10*/; t++)				
						printf("%lf ", hmm.place_gaussian_sigma[i][j][t]);
					printf("\n");
				}

			printf("\nGaussians SIGMA INV SNP %ld:\n", t);
			for(j = 0; j < num_mix_gaussians_sigma; j++)
				for(i = 0; i < num_gaussians; i++)
				{
					for(t = 0; t < seq_len/*10*/; t++)								
						printf("%lf ", hmm.place_gaussian_sigma_inv[i][j][t]);
					printf("\n");
				}
#endif

		}
	}

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached 12\n");
	fclose(debug_f);
#endif
	
	data.seq_len = seq_len; // this does in the data
	data.y_type = CONTINUOUS;
	data.miss_data = miss_data;
	data.dont_miss_prob = 0.5; // not important (only for simulations)

	hmm.seq_len = seq_len;
	hmm.miss_data = miss_data; // double .. 
	hmm.x_dim = x_dim;
	hmm.y_dim = y_dim;
	hmm.y_type = CONTINUOUS;
	hmm.cum_flag = 1;
	hmm.log_flag = 1;
	hmm.place_flag = place_flag; // currently no placing 
	hmm.special_models_flag = special_models_flag;
	hmm.snp_specific_flag = snp_specific_flag;

//	printf("BEFORE INIT MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);




	// Here init the model to get everything in the right place !!! (use_bound is set to 0 since we dont need them in Viterbi)
	InitilizeModel(&hmm, x_dim, y_dim, CONTINUOUS, 0, place_flag, special_models_flag, miss_data, 0);  // fold change flag is set to zero here,,
																		   // since we do not do learning anyway, and don't want to 'ruin' the mus
//	printf("place %d special %d dim2 is %d seq_len %d\n", place_flag, special_models_flag, hmm.x_dim2, hmm.seq_len);	

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached 34\n");
	fclose(debug_f);
#endif

	

	// Now write again the parameters !!!!! why ?? 
	in = mxGetPr(prhs[4]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[4]), mxGetN(prhs[4]));     // Get the x dim
  	for(i = 0; i < x_dim; i++)
		hmm.PI[i] = in[i];   // Copy the initial vector

	in = mxGetPr(prhs[5]);   // Get the transition matrix
	x_dim = mxGetM(prhs[5]);     // Get the x dim
	L = mxGetN(prhs[5]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix



	in = mxGetPr(prhs[6]);   // Get the mixture matrix
	L = mxGetM(prhs[6]);
	y_dim = mxGetN(prhs[6]);     // Get the y dim
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
		if(L != 2*x_dim-1)
		{
			printf("dim: %d 2*dim-1: %d L %d\n", x_dim, 2*x_dim-1, L);
			mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
		}
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.N[i][j] = in[j*(2*x_dim-1)+i];   // Copy the mixture matrix

		in = mxGetPr(prhs[7]);   // Get the mean matrix
	  	for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*(2*x_dim-1)+i];   // Copy the mean matrix

		in = mxGetPr(prhs[8]);   // Get the std matrix
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = MAX(MIN_SIGMA, in[j*(2*x_dim-1)+i]);   // Copy the st.d. matrix

//			for(i = 0; i < 2*x_dim-1; i++)
//			printf("MU %lf SIGMA %lf\n", hmm.MU[i][0], hmm.SIGMA[i][0]);
	
	}
	else
	{
		if(L != x_dim)
		{
			printf("WTFFFFF\n");
			mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
		}
  		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix

		in = mxGetPr(prhs[7]);   // Get the mean matrix
	  	for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

		in = mxGetPr(prhs[8]);   // Get the std matrix
  		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = MAX(MIN_SIGMA, in[j*x_dim+i]);   // Copy the std matrix
	}


	/*****
		in = mxGetPr(prhs[6]);   // Get the mixture matrix
	L = mxGetM(prhs[6]);
	y_dim = mxGetN(prhs[6]);     // Get the y dim
	if(L != x_dim)
		mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[7]);   // Get the mean matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[8]);   // Get the std matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.SIGMA[i][j] = in[j*x_dim+i];   // Copy the std matrix
*****/
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached 56\n");
	fclose(debug_f);
#endif


//	printf("AFTER INIT MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);	
	ComputeModelAuxillaryParameters(&hmm); // Compute the cums ...

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached 78\n");
	fclose(debug_f);
#endif


//	printf("AFTER AUXIL MU12, Sigma12 : %lf %lf | %lf %lf \n", hmm.MU[0][0], hmm.MU[1][0], hmm.SIGMA[0][0], hmm.SIGMA[1][0]);
	// allocate memory


	long max_x_states = x_dim;
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
		max_x_states = 36;	
//		printf("HMM x_dim %ld x_dim2 %ld y_dim %ld y_type %ld\n", hmm.x_dim, hmm.x_dim2, hmm.y_dim, hmm.y_type); 
	/****************************************************************************/
		for(i_geno = 0; i_geno < (hmm.x_dim2*hmm.x_dim2); i_geno++)		// state at time t-1
			for(i1 = 0; i1 < hmm.x_dim; i1++)		// state at time t-1
				for(i2 = 0; i2 < hmm.x_dim; i2++)		// state at time t-1
				{
					// Seperate according to whether the source or destination are zero ..
					total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)] ; 

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Allocate Memory total_i_index %ld\n", total_i_index);
	fclose(debug_f);
#endif

					
					if( (total_i_index <= 4) || (snp_specific_flag && (total_i_index < MAX_2D_X_VALS)) )
					{
						y_cond_tabs[total_i_index] = new double[seq_len];
						y_cond_tabsB[total_i_index] = new double[seq_len];
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, " Pointers Are y_cond_tabs %ld y_cond_tabsB %ld \n", y_cond_tabs[total_i_index], y_cond_tabsB[total_i_index]);
	fclose(debug_f);
#endif

						if( (y_cond_tabs[total_i_index] == NULL) || (y_cond_tabsB[total_i_index] == NULL) )
						{
							printf("Error: failed to allocate hmm.y_cond_tabs\n");
							exit(99);
						}
					}
					alpha[total_i_index] = new double[seq_len];
					beta[total_i_index] = new double[seq_len];
					gamma[total_i_index] = new double[seq_len];
					if(hmm.y_type == CONTINUOUS)
						for(j = 0; j < hmm.y_dim; j++)
						{
							phi[total_i_index][j] = new double[seq_len];
							if( (total_i_index <= 4) || (snp_specific_flag && (total_i_index < MAX_2D_X_VALS)) )
							{
								y_mu_square_exp_vecs[total_i_index][j] = new double[seq_len];
								y_mu_square_exp_vecsB[total_i_index][j] = new double[seq_len];
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, " Pointers Are y_mu_square_exp_vecs %ld y_mu_square_exp_vecsB %ld \n", 
		y_mu_square_exp_vecs[total_i_index][j], y_mu_square_exp_vecsB[total_i_index][j]);
	fclose(debug_f);
#endif

								if( (y_mu_square_exp_vecs[total_i_index][j]  == NULL) || (y_mu_square_exp_vecsB[total_i_index][j] == NULL) )
								{
									printf("Error: failed to allocate hmm.y_cond_tabs\n");
									exit(99);
								}
							}
						}


				}
	/****************************************************************************/
	}
	else  // here special_models_flag != SNPS_ONE_SAMPLE
		for(i = 0; i < hmm.x_dim; i++)	
		{
			// Seperate according to whether the source or destination are zero ..
			y_cond_tabs[i] = new double[seq_len];
			y_cond_tabsB[i] = new double[seq_len];
			alpha[i] = new double[seq_len];
			beta[i] = new double[seq_len];
			gamma[i] = new double[seq_len];
			if(hmm.y_type == CONTINUOUS)
				for(j = 0; j < hmm.y_dim; j++)
				{
					phi[i][j] = new double[seq_len];	
					y_mu_square_exp_vecs[i][j] = new double[seq_len];
					y_mu_square_exp_vecsB[i][j] = new double[seq_len];
				}
		}

	opt_x_vec_out = new long[seq_len];
	scale = new double[seq_len];
	scale_cum = new double[seq_len];
	scale_exp = new double[seq_len];


	long VitInd=2;

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Before VIT\n");
	end_time = clock();
	printf("Beginning SHIT (sec.) : %lf\n", double(end_time-start_time) / double(CLOCKS_PER_SEC));
	start_time = clock();
	printf("THE SPECIAL MODELS FLAG: %ld \n", 	special_models_flag);
	printf("THE PLACE FLAG: %ld\n", place_flag);
	fprintf(debug_f, "POINTER y_cond_tabs 0 NOW %ld\n", y_cond_tabs[0]);
	fclose(debug_f);
#endif
	

//	return; 

//#undef DO_STUFF
#define DO_STUFF
#ifdef DO_STUFF

	// Call many c functions to do the job for you ..
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
#ifdef DEBUG
		start_time = clock();
#endif
		
		long BcondInd=Compute_y_cond_tabsSNPs(&hmm, &data, y_cond_tabs, y_cond_tabsB, 
			y_mu_square_exp_vecs, y_mu_square_exp_vecsB);  // OK
#ifdef DEBUG
		end_time = clock();
#endif	

	}
	else
		Compute_y_cond_tabs(&hmm, &data, y_cond_tabs, y_mu_square_exp_vecs);  // OK, old models

#ifdef DEBUG
	end_time = clock();
	printf("B Cond TAB (sec.) : %lf\n", double(end_time-start_time) / double(CLOCKS_PER_SEC));
#endif


///////////////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************** All Phases *******/


#ifdef DEBUG
	// Print some of them
	long A_copy=0, B_copy=1, tt=4;
	printf("On Position Four (4): A_copy = 0 B_copy = 1: \n");
	printf("mu_A %lf  mu_B %lf sigma_A %lf sigma_B %lf cov_AB %lf sigma_inv_A %lf sigma_inv_B %lf sigma_inv_AB %lf\n", 
		hmm.place_gaussian_mu[A_copy+5*B_copy][0][tt], hmm.place_gaussian_mu[A_copy+5*B_copy][1][tt], 
		hmm.place_gaussian_sigma[A_copy+5*B_copy][0][tt], hmm.place_gaussian_sigma[A_copy+5*B_copy][1][tt], 
		hmm.place_gaussian_sigma[A_copy+5*B_copy][2][tt], 
		hmm.place_gaussian_sigma_inv[A_copy+5*B_copy][0][tt], hmm.place_gaussian_sigma_inv[A_copy+5*B_copy][1][tt], 
		hmm.place_gaussian_sigma_inv[A_copy+5*B_copy][2][tt]); 


	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Y EXP VEC (A Copy ...)\n");
	for(A_copy = 0; A_copy <= 4; A_copy++)		
		for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
		{
			for(t=0; t < /*seq_len*/ 20; t++)
				fprintf(debug_f, "%lf, ", y_mu_square_exp_vecs[A_copy+5*B_copy][0][t]);
			fprintf(debug_f, "\n");
		}

	fprintf(debug_f, "Y COND (A Copy ...)\n");
	for(A_copy = 0; A_copy <= 4; A_copy++)		
		for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
		{
			fprintf(debug_f, "A Copy %ld B Copy %ld | ", A_copy, B_copy);
			for(t=0; t < /*seq_len*/ 20; t++)
				fprintf(debug_f, "%lf, ", y_cond_tabs[A_copy+5*B_copy][t]);
			fprintf(debug_f, "\n");
		}



	fprintf(debug_f, "Data (A Copy ...)\n");
	for(A_copy = 0; A_copy <= 4; A_copy++)		
		for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
		{
			fprintf(debug_f, "A Copy %ld B Copy %ld | ", A_copy, B_copy);
			for(t=0; t < /*seq_len*/ 20; t++)
				fprintf(debug_f, "%lf, ", data.y_vec[t]  );
			fprintf(debug_f, "\n");
		}
	fprintf(debug_f, "MU (A Copy ...)\n");
	for(A_copy = 0; A_copy <= 4; A_copy++)		
		for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
		{
			fprintf(debug_f, "A Copy %ld B Copy %ld | ", A_copy, B_copy);
			for(t=0; t < /*seq_len*/ 20; t++)
				fprintf(debug_f, "%lf, ", hmm.place_gaussian_mu[A_copy+5*B_copy][0][t] );
			fprintf(debug_f, "\n");
		}



	fprintf(debug_f, "SQR DIFF (A Copy ...)\n");
	for(A_copy = 0; A_copy <= 4; A_copy++)		
		for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
		{
			fprintf(debug_f, "A Copy %ld B Copy %ld | ", A_copy, B_copy);
			for(t=0; t < /*seq_len*/ 20; t++)
				fprintf(debug_f, "%lf, ", (data.y_vec[t]-hmm.place_gaussian_mu[A_copy+5*B_copy][0][t])* 
				(data.y_vec[t]-hmm.place_gaussian_mu[A_copy+5*B_copy][0][t])   );
			fprintf(debug_f, "\n");
		}
	fprintf(debug_f, "SQR DIFF (B Copy ...)\n");
	for(A_copy = 0; A_copy <= 4; A_copy++)		
		for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
		{
			fprintf(debug_f, "A Copy %ld B Copy %ld | ", A_copy, B_copy);
			for(t=0; t < /*seq_len*/ 20; t++)
				fprintf(debug_f, "%lf, ", (data.y_vecB[t]-hmm.place_gaussian_mu[A_copy+5*B_copy][1][t])* 
				(data.y_vecB[t]-hmm.place_gaussian_mu[A_copy+5*B_copy][1][t])   );
			fprintf(debug_f, "\n");
		}


	fprintf(debug_f, "SQR DIFF (SSQ Copy ...)\n");
	for(A_copy = 0; A_copy <= 4; A_copy++)		
		for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
		{
			fprintf(debug_f, "A Copy %ld B Copy %ld | ", A_copy, B_copy);
			for(t=0; t < /*seq_len*/ 20; t++)
				fprintf(debug_f, "%lf, ", (data.y_vec[t]-hmm.place_gaussian_mu[A_copy+5*B_copy][0][t])* 
				(data.y_vec[t]-hmm.place_gaussian_mu[A_copy+5*B_copy][0][t]) + 
				(data.y_vecB[t]-hmm.place_gaussian_mu[A_copy+5*B_copy][1][t])* 
				(data.y_vecB[t]-hmm.place_gaussian_mu[A_copy+5*B_copy][1][t]) );
			fprintf(debug_f, "\n");
		}



/***
	printf("Y CONDB (B Copy ...)\n");
	for(i=0; i <= 4; i++)
	{
		for(t=0; t < 10; t++)
			printf("%lf, ", y_cond_tabsB[i][t]);
		printf("\n");
	}
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "BCOND IS %d\n", VitInd);
***/
	fclose(debug_f);
#endif


	double x_loglike, x_place_loglike, y_loglike1, y_loglike2;
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
#ifdef DEBUG
		debug_f = fopen(debug_p, "a");
		fprintf(debug_f, "BEFORE Forward SNPs - The INPUTS Are: \n");


		// Print some of the first data 			
		fprintf(debug_f, "data y :\n", t);
		for(t=0; t < seq_len/*10*/; t++)
			fprintf(debug_f, "%lf, ", data.y_vec[t]);
		fprintf(debug_f, "\n");
			
		fprintf(debug_f, "data yB :\n", t);
		for(t=0; t < seq_len/*10*/; t++)
			fprintf(debug_f, "%lf, ", data.y_vecB[t]);
		fprintf(debug_f, "\n");

		fprintf(debug_f, "HMM Params\n"); 
		PrintModel2(&hmm);
		
		fprintf(debug_f, "Emission Matrix: \n");
		for(i = 0; i < hmm.x_dim; i++)
		{
			for(j = 0; j < hmm.x_dim; j++)
				fprintf(debug_f, "%lf, ", hmm.M[i][j]);
			fprintf(debug_f, "\n");
		}

		fprintf(debug_f, "x_dim %ld x_dim2 %ld seq_len %ld gauss_dim %ld\n", hmm.x_dim, hmm.x_dim2, data.seq_len, hmm.gauss_dim);


		fprintf(debug_f, "PLACE M matrix: \n");

		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
			{
				for(t=0; t < seq_len/*10*/; t++)
					fprintf(debug_f, "%lf ", hmm.place_M[i][j][t]); 
				fprintf(debug_f, "\n");
			}


		
		// Print some of the first y_conds 			
		fprintf(debug_f, "y_cond_tabs :\n", t);
		for(A_copy = 0; A_copy <= 4; A_copy++)		
			for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
			{
				fprintf(debug_f, "A Copy %ld B Copy %ld | ", A_copy, B_copy);
				for(t=0; t < /*seq_len*/ 20; t++)
					fprintf(debug_f, "%lf, ", y_cond_tabs[A_copy+5*B_copy][t]);
				fprintf(debug_f, "\n");
			}
		fclose(debug_f);
		start_time = clock();
#endif

//		return; // early termination - after forward


		forwardSNPs(&hmm, &data, y_cond_tabs, y_cond_tabsB, alpha, scale, scale_cum, scale_exp, &cur_score);  // PROBLEM FIRST ALPHA !!! 

	
#ifdef DEBUG
		end_time = clock();
		printf("Forward SNPs (sec.) : %lf\n", double(end_time-start_time) / double(CLOCKS_PER_SEC));
#endif

//		return; // early termination - after forward


#ifdef DEBUG	
		debug_f = fopen(debug_p, "a");
		fprintf(debug_f, "After Forward SNPs\n");

		// Print some of the first alpha's 			
		fprintf(debug_f, "Alpha values :\n", t);
		for(i = 0; i < 36; i++)
		{
			for(t = 0; t < /*seq_len*/10; t++)
				fprintf(debug_f, "%lf ", alpha[i][t]);
			printf("\n\n");
		}

		
		fclose(debug_f);
		start_time = clock();
#endif
		backwardSNPs(&hmm, &data, scale_exp, y_cond_tabs, y_cond_tabsB, beta);		

//		return; // early termination - after forward
#ifdef DEBUG	

		end_time = clock();
		printf("Backward SNPs (sec.) : %lf\n", double(end_time-start_time) / double(CLOCKS_PER_SEC));
/***	printf("OLD BETA\n");
		for(i=0; i < 10; i++)
		{
			for(j=0; j < 36; j++)
				printf("%lf ", beta[i][j]);
			printf("\n");
		} ***/
#endif



#ifdef DEBUG
		debug_f = fopen(debug_p, "a");
		fprintf(debug_f, "Backwards SNPs\n");

		// Print some of the first beta's 			
		fprintf(debug_f, "\nBeta values :\n", t);
		for(i = 0; i < 36; i++)
		{
			for(t = 0; t < /*seq_len*/10; t++)
				fprintf(debug_f, "%lf ", beta[i][t]);
			printf("\n\n");
		}

		fclose(debug_f);
		start_time = clock();
#endif
		ComputeMarginalGammaSNPs(&hmm, &data, y_mu_square_exp_vecs, y_mu_square_exp_vecsB, 
								 alpha, beta, scale_exp, y_cond_tabs, y_cond_tabsB,
								 gamma, phi);  // phi is not used here ...
#ifdef DEBUG
		end_time = clock();
		printf("MarginalGamma SNPs (sec.) : %lf\n", double(end_time-start_time) / double(CLOCKS_PER_SEC));
		debug_f = fopen(debug_p, "a");
		fprintf(debug_f, "GammaSNPs \n");

				// Print some of the first beta's 			
		fprintf(debug_f, "\nGamma values :\n", t);
		for(i = 0; i < 36; i++)
		{
			for(t = 0; t < /*seq_len*/10; t++)
				fprintf(debug_f, "%lf ", gamma[i][t]);
			fprintf(debug_f, "\n");
		}


		fclose(debug_f);
		start_time = clock();
#endif

//		return; // early termination - after forward


		VitInd = ViterbiSNPs( &hmm, &data, y_cond_tabs, y_cond_tabsB, opt_x_vec_out); // , &x_loglike, &x_place_loglike, &y_loglike1);

//		return; // early termination - after forward

#ifdef DEBUG
		end_time = clock();
		printf("Viterbi SNPs (sec.) : %lf\n", double(end_time-start_time) / double(CLOCKS_PER_SEC));
		printf("Viterbi Path\n");
//		for(i=0; i < 50; i++)
//			printf("%ld ", opt_x_vec_out[i]); 
//		printf("\n");
		debug_f = fopen(debug_p, "a");
		fprintf(debug_f, "VitSNPs %ld \n", VitInd);
		fclose(debug_f);
#endif

#ifdef DEBUG
		debug_f = fopen(debug_p, "a");
		fprintf(debug_f, "VitSNPs %ld \n", VitInd);
		fclose(debug_f);
		start_time = clock();
#endif



		// Add to check the likelihood ...
/***
		double cur_score, x_loglike, *place_x_loglike, y_loglike;
		data.x_vec = new long[seq_len];
		place_x_loglike = new double[seq_len];
		for(i=0; i<seq_len; i++)
			data.x_vec[i] = opt_x_vec_out[i];
		cur_score = ComputeJointHMMLogLikelihoodSNPs(&hmm, &data, &x_loglike, place_x_loglike, &y_loglike);
***/

/*****
		printf("V LogLike parts: x_loglike %lf place_x_loglike %lf y_loglike %lf ALL %lf\n", 
			x_loglike, place_x_loglike[0], y_loglike, cur_score);
		for(i=0; i<10; i++)
			printf("I %ld PLACE LL %lf\n", i, place_x_loglike[i]);

/********
		double pi_vec[MAX_X_VALS];
		MatrixStationaryVec(hmm.M, hmm.x_dim, pi_vec);

		// print both of them: 
		printf("The Markov Matrix: \n");
		for(i=0;i<hmm.x_dim;i++)
		{
			for(j=0;j<hmm.x_dim;j++)
				printf("%lf ", hmm.M[i][j]);
			printf("\n");
		}
		printf("The Stationary Vector: \n");
		for(i=0;i<hmm.x_dim;i++)
			printf("%lf ", pi_vec[i]);
		printf("\n");
********/
	}
	else
	{
#ifdef DEBUG
		debug_f = fopen(debug_p, "a");
		fprintf(debug_f, "BEFORE Forward \n");
		fclose(debug_f);
#endif
		forward(&hmm, &data, y_cond_tabs, alpha, scale, scale_cum, scale_exp, &cur_score);  // PROBLEM FIRST ALPHA !!! 
#ifdef DEBUG
		debug_f = fopen(debug_p, "a");
		fprintf(debug_f, "Forward \n");
		fclose(debug_f);
#endif
		backward(&hmm, &data, scale_exp, y_cond_tabs, beta);		
		ComputeMarginalGamma(&hmm, &data, y_mu_square_exp_vecs, alpha, beta, scale_exp, y_cond_tabs,
							gamma, phi);
		VitInd = Viterbi( &hmm, &data, y_cond_tabs, y_cond_tabs, opt_x_vec_out, &x_loglike, &x_place_loglike, &y_loglike1);
	}

	// Now copy Viterbi output vector  : 
	plhs[0] = mxCreateDoubleMatrix(1, seq_len, mxREAL);  // The Viterbi vector 
	out = mxGetPr(plhs[0]); // output Viterbi
	for(i=0; i<seq_len; i++) 
		out[i] = opt_x_vec_out[i];

	/****
	printf("The X_dim %ld The Viterbi output:\n", hmm.x_dim);
	for(i=0; i<seq_len; i++) 
		printf("%ld ", opt_x_vec_out[i]);
	printf("\n");
****/

//	return;

	// Now copy gamma probs output vector
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
		long TotalDim = hmm.x_dim*hmm.x_dim*4; // contains both genotypes and both copynumbers
		plhs[1] = mxCreateDoubleMatrix(TotalDim, seq_len, mxREAL);  // The gamma vector 
		out = mxGetPr(plhs[1]); // output gamma
		for(i_geno = 0; i_geno < (hmm.x_dim2*hmm.x_dim2); i_geno++)		// state at time t-1
			for(i1 = 0; i1 < hmm.x_dim; i1++)		// state at time t-1
				for(i2 = 0; i2 < hmm.x_dim; i2++)		// state at time t-1
				{
					j=BIT(i_geno,0) + (i1 << 1) + 2*hmm.x_dim*BIT(i_geno,1) + 4*hmm.x_dim*i2; // here it is compressed
					// Seperate according to whether the source or destination are zero ..
					total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)] ; 
					for(t=0; t<seq_len; t++) 
						out[j+t*TotalDim] = gamma[total_i_index][t]; // should be gamma but we want to check and debug .. .
				}

		// WRRRONGGG: FOR DEBUGGING:
		/***
		for(i = 0; i < 5; i++)
			for(t=0; t<seq_len; t++) 
			{
				out[i+t*TotalDim] = y_cond_tabs[i][t]; // should be Prob Y_A Given MU_A  .. .
				out[i+10+t*TotalDim] = y_cond_tabsB[i][t]; // should be Prob Y_B Given MU_B  .. .
			}
			for(i=1; i < 1000; i += 10)
				printf("I %ld  Exp(-I) %lf\n", i, 1.0 / exp(-i)); // *11111111111111111111.111);
				/***/


	}
	else
	{
		plhs[1] = mxCreateDoubleMatrix(x_dim, seq_len, mxREAL);  // The gamma vector 
		out = mxGetPr(plhs[1]); // output gamma
		for(i = 0;i < x_dim; i++)
			for(t=0; t<seq_len; t++) 
				out[i+t*x_dim] = gamma[i][t]; // gamma[i][t];
	}

#endif // DO_STUFF

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached OUTPUT Now X_DIM %ld Y_DIM %ld snp_specific_flag %ld\n", hmm.x_dim, hmm.y_dim, snp_specific_flag);
	fprintf(debug_f, "POINTER y_cond_tabs 0 NNNOW %ld\n", y_cond_tabs[0]);
	fclose(debug_f);
#endif


//	return; // here only printed something

/************************************************************************************** End All Phases ********/
///////////////////////////////////////////////////////////////////////////////////////////////////////////

	/******* FREE MEMORY ****/
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
		max_x_states = 36;	
	/****************************************************************************/
		for(i_geno = 0; i_geno < (hmm.x_dim2*hmm.x_dim2); i_geno++)		// state at time t-1
			for(i1 = 0; i1 < hmm.x_dim; i1++)		// state at time t-1
				for(i2 = 0; i2 < hmm.x_dim; i2++)		// state at time t-1
				{

					// Seperate according to whether the source or destination are zero ..
					total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)] ; 
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	if(total_i_index == 0)
		fprintf(debug_f, "POINTER y_cond_tabs 0 NNNOW %ld\n", y_cond_tabs[0]);
	fclose(debug_f);
#endif
					

					if( (total_i_index <= 4) || (snp_specific_flag && (total_i_index < MAX_2D_X_VALS)) )
					{
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Free Memory Cond Tabs total_i_index %ld\n", total_i_index);
	fclose(debug_f);
#endif

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, " Pointers Are y_cond_tabs %ld y_cond_tabsB %ld \n", y_cond_tabs[total_i_index], y_cond_tabsB[total_i_index]);
	fclose(debug_f);
#endif

						delete y_cond_tabs[total_i_index];
						delete y_cond_tabsB[total_i_index];
					}

					delete alpha[total_i_index];
					delete beta[total_i_index];
					delete gamma[total_i_index];
/***/
					if(hmm.y_type == CONTINUOUS)
						for(j = 0; j < hmm.y_dim; j++)
						{
							delete phi[total_i_index][j];	
							if( (total_i_index <= 4) || (snp_specific_flag && (total_i_index < MAX_2D_X_VALS)) )
							{
								delete y_mu_square_exp_vecs[total_i_index][j];
								delete y_mu_square_exp_vecsB[total_i_index][j];
							}
						}
						/***/
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Free Memory exp vec total_i_index %ld\n", total_i_index);
	fclose(debug_f);
#endif

				}
//		return; 	

	/****************************************************************************/
	}
	else
		for(i = 0; i < hmm.x_dim; i++)
		{
			delete y_cond_tabs[i];
			delete y_cond_tabsB[i];
			delete alpha[i];
			delete beta[i];
			delete gamma[i];
			for(j = 0; j < hmm.y_dim; j++)
			{
				delete y_mu_square_exp_vecs[i][j];
				delete y_mu_square_exp_vecsB[i][j];
				delete phi[i][j];
			}
		}

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Middle Free Memeory\n");
	fclose(debug_f);
#endif

//	return; 

	// delete memory 
	delete data.y_vec;
	delete data.y_vecB;
	delete data.loc_vec;
	delete data.loc_diff_vec;	
	delete opt_x_vec_out;
	delete scale;
	delete scale_cum; 
	delete scale_exp;
	if(place_flag)
	{	
		if(special_models_flag)
			for(i = 0; i < hmm.x_dim2; i++)		
				for(j = 0; j < hmm.x_dim2; j++)
				{
					delete hmm.place_M[i][j];
					delete hmm.place_M_cum[i][j];
				}
		if(snp_specific_flag || (~special_models_flag))
			for(i = 0; i < num_gaussians; i++)		
			{
				for(j = 0; j < num_mix_gaussians_mu; j++)
					delete hmm.place_gaussian_mu[i][j];
				for(j = 0; j < num_mix_gaussians_sigma; j++)
				{
					delete hmm.place_gaussian_sigma[i][j];
					delete hmm.place_gaussian_sigma_inv[i][j];
				}
			}
	}

//	return;

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Finished Free Memeory\n");
	fclose(debug_f);
	end_time = clock();
	printf("END STUFF SHIT (sec.) : %lf\n", double(end_time-start_time) / double(CLOCKS_PER_SEC));

#endif
	

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






