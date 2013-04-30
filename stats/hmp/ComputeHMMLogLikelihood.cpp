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
  

	// Open again and start reading
	char *debug_p = "DEBUG.txt";
	FILE *debug_f;
	debug_f = fopen(debug_p, "w");
	fprintf(debug_f, "Starting\n");
	fclose(debug_f);

	long nRows, nCols, i, j, K, t, max_iters, num_starts;
	long i1,i2,j1,j2,i_geno,j_geno,total_i_index;

	long use_x_flag, seq_len, L; 
	long x_dim, x_dim2, y_dim, place_flag, miss_data, special_models_flag;
	double tolerance;
  
	double *in, *out;

	hmm_model hmm;  // containts the model
	hmm_data data; // contains the data

	double *y_cond_tabs[MAX_2D_X_VALS]; // helpful conditional tables
	double *y_cond_tabsB[MAX_2D_X_VALS]; // helpful conditional tables
	double *alpha[36];
	double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS];
	double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS];

	double *scale, *scale_cum, *scale_exp;
	double cur_score;

	 /* Check for proper number of arguments. Should be Fourteen 14 */
	if(nrhs != 16) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. Hidden State vector (X)
	2. Expression vector
	3. Second Expression vector
	4. Flag saying if to use X
	5. Chromosomal location vector
	6. missing data flag 
	7. PI initial vector
	8. M transition matrix
	9. N mixture matrix 
	10. MU mean matrix
	11. SIGMA standard error matrix
	12. place flag
	13. special model flag (NEW)	
	14. place means
	15. place sigmas
	16. place M transition matrices

 ---   Note : Dimensions are already IN the model !!! ---
	*/
	in = mxGetPr(prhs[3]);   // Get the use_x flag
	use_x_flag = long(in[0]);  
	

	if(use_x_flag)
	{
		in = mxGetPr(prhs[0]);   // Get the expression vector
		seq_len = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the expression vector length
		data.x_vec = new long[seq_len];
		for(t = 0; t < seq_len; t++)   // copy expression vector
			data.x_vec[t] = long(in[t]);  
	}


	in = mxGetPr(prhs[1]);   // Get the expression vector
	seq_len = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));     // Get the expression vector length
	data.y_vec = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy expression vector
		data.y_vec[t] = in[t];  

	in = mxGetPr(prhs[2]);   // Get the additional expression vector if needed
	data.y_vecB = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy expression vector
		data.y_vecB[t] = in[t];  

	in = mxGetPr(prhs[4]);   // Get the location vector
	L  = MAX(mxGetM(prhs[4]), mxGetN(prhs[4]));    // Get the number of sequences (genes)
	if(L != seq_len)
		mexErrMsgTxt("Usage: Exp. and Loc. vector length must be the same");	
	data.loc_vec = new long[seq_len];
	for(t = 0; t < seq_len; t++)   // copy location vector
		data.loc_vec[t] = long(in[t]/NUM_BASEPAIRS_PER_HMM_TRANS);  
	data.loc_diff_vec = new long[seq_len];
	data.loc_diff_vec[0] = 0; 
	for(t = 1; t < seq_len; t++)   // compute location diff vector
		data.loc_diff_vec[t] = MIN(MAX(data.loc_vec[t]-data.loc_vec[t-1]-1, 0), MAX_POWER_OF_TRANS_MATRIX-2);  


	in = mxGetPr(prhs[5]);   // Get the missing data flag (use locations)
	miss_data = long(in[0]);  
	in = mxGetPr(prhs[6]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[6]), mxGetN(prhs[6]));     // Get the x dim

  	for(i = 0; i < x_dim; i++)
		hmm.PI[i] = in[i];   // Copy the initial vector


	in = mxGetPr(prhs[7]);   // Get the transition matrix
	x_dim = mxGetM(prhs[7]);     // Get the x dim
	L = mxGetN(prhs[7]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix


	in = mxGetPr(prhs[11]);   // Get the place flag
	place_flag = long(in[0]);
	in = mxGetPr(prhs[12]);   // Get the special models flag
	special_models_flag = long(in[0]);


	in = mxGetPr(prhs[8]);   // Get the mixture matrix
	L = mxGetM(prhs[8]);
	y_dim = mxGetN(prhs[8]);     // Get the y dim
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

		in = mxGetPr(prhs[9]);   // Get the mean matrix
	  	for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*(2*x_dim-1)+i];   // Copy the mean matrix

		in = mxGetPr(prhs[10]);   // Get the std matrix
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = MAX(MIN_SIGMA, in[j*(2*x_dim-1)+i]);   // Copy the st.d. matrix
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

		in = mxGetPr(prhs[9]);   // Get the mean matrix
	  	for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

		in = mxGetPr(prhs[10]);   // Get the std matrix
  		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = MAX(MIN_SIGMA, in[j*x_dim+i]);   // Copy the std matrix. Avoid too little std.s.
	}


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
						hmm.place_M[i][j][t] = in[t*x_dim2*x_dim2+j*x_dim2+i]; // in[i+x_dim2*j+x_dim2*x_dim2*t]; // Copy the place means

		}
		else
		{
			// allocate memory 
			for(i = 0; i < x_dim; i++)		
				for(j = 0; j < y_dim; j++)
				{
					hmm.place_gaussian_mu[i][j] = new double[seq_len];
					hmm.place_gaussian_sigma[i][j] = new double[seq_len];
				}

			// Now copy
			in = mxGetPr(prhs[13]);   // Get the place means
			/////long arr_size = mxGetM(prhs[10]) * mxGetN(prhs[10]);
			for(i = 0; i < x_dim; i++)
				for(j = 0; j < y_dim; j++)
					for(t = 0; t < seq_len; t++)
						hmm.place_gaussian_mu[i][j][t] = in[t+seq_len*j+seq_len*y_dim*i]; // Copy the place means
			// Now copy
			in = mxGetPr(prhs[14]);   // Get the place sigmas
			for(i = 0; i < x_dim; i++)	
				for(j = 0; j < y_dim; j++)
					for(t = 0; t < seq_len; t++)
						hmm.place_gaussian_sigma[i][j][t] = in[t+seq_len*j+seq_len*y_dim*i]; // Copy the place sigmas
		}
	}

	
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


	// Here init the model to get everything in the right place !!! (use_bound is set to 0 since we dont need them in Viterbi)
	InitilizeModel(&hmm, x_dim, y_dim, CONTINUOUS, 0, place_flag, special_models_flag, miss_data, 0);  // fold change flag is set to zero here,,
																		   // since we do not do learning anyway, and don't want to 'ruin' the mus
	

	// Now write again the parameters !!!!! why ?? 
	in = mxGetPr(prhs[6]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[6]), mxGetN(prhs[6]));     // Get the x dim
  	for(i = 0; i < x_dim; i++)
		hmm.PI[i] = in[i];   // Copy the initial vector

	in = mxGetPr(prhs[7]);   // Get the transition matrix
	x_dim = mxGetM(prhs[7]);     // Get the x dim
	L = mxGetN(prhs[7]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix
	in = mxGetPr(prhs[8]);   // Get the mixture matrix
	L = mxGetM(prhs[8]);
	y_dim = mxGetN(prhs[8]);     // Get the y dim
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

		in = mxGetPr(prhs[9]);   // Get the mean matrix
	  	for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*(2*x_dim-1)+i];   // Copy the mean matrix

		in = mxGetPr(prhs[10]);   // Get the std matrix
  		for(i = 0; i < 2*x_dim-1; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = MAX(MIN_SIGMA, in[j*(2*x_dim-1)+i]);   // Copy the st.d. matrix
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

		in = mxGetPr(prhs[9]);   // Get the mean matrix
	  	for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.MU[i][j] = in[j*x_dim+i];   // Copy the mean matrix

		in = mxGetPr(prhs[10]);   // Get the std matrix
  		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				hmm.SIGMA[i][j] = MAX(MIN_SIGMA, in[j*x_dim+i]);   // Copy the std matrix
	}

	ComputeModelAuxillaryParameters(&hmm); // Compute the cums ...


	// allocate memory
	long max_x_states = x_dim;
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
		max_x_states = 36;	
	//	printf("HMM x_dim %ld x_dim2 %ld y_dim %ld y_type %ld\n", hmm.x_dim, hmm.x_dim2, hmm.y_dim, hmm.y_type); 
	/****************************************************************************/
		for(i_geno = 0; i_geno < (hmm.x_dim2*hmm.x_dim2); i_geno++)		// state at time t-1
			for(i1 = 0; i1 < hmm.x_dim; i1++)		// state at time t-1
				for(i2 = 0; i2 < hmm.x_dim; i2++)		// state at time t-1
				{
					// Seperate according to whether the source or destination are zero ..
					total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 
					alpha[total_i_index] = new double[seq_len];
					if(total_i_index <= 4)
					{
						y_cond_tabs[total_i_index] = new double[seq_len];
						y_cond_tabsB[total_i_index] = new double[seq_len];
						if(hmm.y_type == CONTINUOUS)
							for(j = 0; j < hmm.y_dim; j++)
							{
								y_mu_square_exp_vecs[total_i_index][j] = new double[seq_len];
								y_mu_square_exp_vecsB[total_i_index][j] = new double[seq_len];
							}
					}
				}
	/****************************************************************************/
	}
	else
		for(i = 0; i < hmm.x_dim; i++)	
		{
			// Seperate according to whether the source or destination are zero ..
			y_cond_tabs[i] = new double[seq_len];
			y_cond_tabsB[i] = new double[seq_len];
			alpha[i] = new double[seq_len];
			if(hmm.y_type == CONTINUOUS)
				for(j = 0; j < hmm.y_dim; j++)
				{
					y_mu_square_exp_vecs[i][j] = new double[seq_len];
					y_mu_square_exp_vecsB[i][j] = new double[seq_len];
				}
		}


	scale = new double[seq_len];
	scale_cum = new double[seq_len];
	scale_exp = new double[seq_len];


	long VitInd=2;

	double x_loglike, *place_x_loglike, y_loglike;

	// Call many c functions to do the job for you ..
	if(special_models_flag == SNPS_ONE_SAMPLE)
	{
		if(use_x_flag)
		{
			place_x_loglike = new double[seq_len];
			cur_score = ComputeJointHMMLogLikelihoodSNPs(&hmm, &data, &x_loglike, place_x_loglike, &y_loglike);
		}
		else
		{
			long BcondInd=Compute_y_cond_tabsSNPs(&hmm, &data, y_cond_tabs, y_cond_tabsB, y_mu_square_exp_vecs, y_mu_square_exp_vecsB);  // OK
			forwardSNPs(&hmm, &data, y_cond_tabs, y_cond_tabsB, alpha, scale, scale_cum, scale_exp, &cur_score);  // PROBLEM FIRST ALPHA !!! 
			
/***
			printf("First few Alpha's\n");
			for(i=0;i<10;i++)
				printf("%lf, ", alpha[0][i]);
			printf("\n");
***/
		}

//		printf("LogLike parts: x_loglike %lf place_x_loglike %lf y_loglike %lf\n", x_loglike, place_x_loglike[0], y_loglike);
/***
		for(i=0;i<10;i++)
			printf("I %ld PLACE LOG LIKE %lf\n", i, place_x_loglike[i]);
***/
	}
	else
	{
				printf("THIS IS NOT SPECIAL\n");
		if(use_x_flag)
		{
			printf("AM I in the right place?\n");
			cur_score = ComputeJointHMMLogLikelihood(&hmm, &data);
/***
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			double loglike_score = 0.0; // This is what we return
			double log_x = 0.0; 
			double log_y = 0.0; 
			
			long t;
			double pi = 3.1415927; 

			// First get the score for the X's 
			for(t = 0; t < data.seq_len-1; t++)
			{				
				loglike_score += log( hmm.M[data.x_vec[t]][data.x_vec[t+1]] ) - 
					(data.y_vec[t] - hmm.MU[data.x_vec[t]][0]) * (data.y_vec[t] - hmm.MU[data.x_vec[t]][0]) / 
					(2*hmm.SIGMA[data.x_vec[t]][0]*hmm.SIGMA[data.x_vec[t]][0]) - 
					0.5*log(2*pi) - log(hmm.SIGMA[data.x_vec[t]][0]); 
				log_x += log( hmm.M[data.x_vec[t]][data.x_vec[t+1]] );
				log_y += ( -(data.y_vec[t] - hmm.MU[data.x_vec[t]][0]) * (data.y_vec[t] - hmm.MU[data.x_vec[t]][0]) / 
					(2*hmm.SIGMA[data.x_vec[t]][0]*hmm.SIGMA[data.x_vec[t]][0]) - 
					0.5*log(2*pi) - log(hmm.SIGMA[data.x_vec[t]][0]) ); 

				if(t < 10)
					printf("x %ld  y %lf  mu %lf sigma %lf lll_x %lf  lll_y %lf lll_rest_y %lf\n", 
					data.x_vec[t], data.y_vec[t], hmm.MU[data.x_vec[t]][0], hmm.SIGMA[data.x_vec[t]][0], 
					log( hmm.M[data.x_vec[t]][data.x_vec[t+1]] ),
					(data.y_vec[t] - hmm.MU[data.x_vec[t]][0]) * (data.y_vec[t] - hmm.MU[data.x_vec[t]][0]) / 
						(2*hmm.SIGMA[data.x_vec[t]][0]*hmm.SIGMA[data.x_vec[t]][0]), 
											-0.5*log(2*pi) - log(hmm.SIGMA[data.x_vec[t]][0]));

			}
			printf("Output Scores:log_x_score %lf log_y_score %lf\n", log_x / (data.seq_len-1) , log_y / (data.seq_len-1));
***/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		}
		else
		{
			Compute_y_cond_tabs(&hmm, &data, y_cond_tabs, y_mu_square_exp_vecs);  // OK, old models
			forward(&hmm, &data, y_cond_tabs, alpha, scale, scale_cum, scale_exp, &cur_score);  // PROBLEM FIRST ALPHA !!! 
		}
	}

//	printf("Output Scores: scale %lf scale_cum %lf scale_exp %lf cur_score %lf\n", 
//		*scale, *scale_cum, *scale_exp, cur_score);




	// Copy output: Scale exponent: 
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);  // The Scale exponent
	out = mxGetPr(plhs[0]); // output Viterbi
	out[0] = cur_score;


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
					total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 				
					delete alpha[total_i_index];
					if(total_i_index <= 4)
					{
						delete y_cond_tabs[total_i_index];
						delete y_cond_tabsB[total_i_index];
						if(hmm.y_type == CONTINUOUS)
							for(j = 0; j < hmm.y_dim; j++)
							{
								delete y_mu_square_exp_vecs[total_i_index][j];
								delete y_mu_square_exp_vecsB[total_i_index][j];
							}
					}
				}
	/****************************************************************************/
	}
	else
		for(i = 0; i < hmm.x_dim; i++)
		{
			delete y_cond_tabs[i];
			delete y_cond_tabsB[i];
			delete alpha[i];
			for(j = 0; j < hmm.y_dim; j++)
			{
				delete y_mu_square_exp_vecs[i][j];
				delete y_mu_square_exp_vecsB[i][j];
			}
		}

	// delete memory 
	if(use_x_flag)
		delete data.x_vec;
	delete data.y_vec;
	delete data.y_vecB;
	delete data.loc_vec;
	delete data.loc_diff_vec;	
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
		else
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






