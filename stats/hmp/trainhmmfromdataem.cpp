#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
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
#ifdef DEBUG
	char *debug_p = "TrainDEBUG.txt";
	FILE *debug_f;
	debug_f = fopen(debug_p, "w");
	fprintf(debug_f, "Starting\n");
	fclose(debug_f);
#endif
	
	int nRows, nCols, i, j, K, t, max_iters, num_starts;

	long seq_len, L; 
	long x_dim, y_dim, place_flag, miss_data, fold_change_flag;
	double tolerance;
	double *in, *out;
	double best_model_score;

	hmm_model hmm;   // a struct containing the hmm parameters
	hmm_data data;  // a struct containing all the data needed


	 /* Check for proper number of arguments. Should be Seventeen 17 */
	if(nrhs != 18) 
		mexErrMsgTxt("Usage: [dat] [K] 18 parameters");

  /* Parameters should be : 
	1. Expression vector
	2. Expression vector 2
	3. Chromosomal location vector
	4. Missing data flag (says if to use locations ... )
	5. Starts of segments
	6. X dimension required
	7. Y dimension (mixtures) required
	8. place flag (if to use different distribution for each place)
	9. fold change flag (if to bound the mus by one) 
	10. mean vec (currently it is given and not learned)
	11. std vec (currently it is given and not learned)
	12. place_M (The place matrices, should be 2*2 from the HAPMAP) - we don't learn them! 
	13. Special SNP model flag. If this is off we do 'standard' EM. If its on we do the special SNPs model
	14. UpperBounds matrix (bounds on entries of the transition matrix)
	15. bound flag (says if to use upper bounds)
	16. number of EM iterations
	17. Number of EM starting points 
	18. EM tolerance 

	We also need to READ the place_M matrices (they're not learned from the data ...)
	*/

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached 0\n");
	fclose(debug_f);
#endif

	in = mxGetPr(prhs[0]);   // Get the expression vector
	seq_len = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the expression vector length
	data.y_vec = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy expression vector
		data.y_vec[t] = in[t];  
	in = mxGetPr(prhs[1]);   // Get the expression vector number 2
	seq_len = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));     // Get the expression vector2 length
	data.y_vecB = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy expression vector number 2
		data.y_vecB[t] = in[t];  
	in = mxGetPr(prhs[2]);   // Get the location vector
	L  = MAX(mxGetM(prhs[2]), mxGetN(prhs[2]));    // Get the location vector length
	if(L != seq_len)
		mexErrMsgTxt("Usage: Exp. and Loc. vector length must be the same");
	data.loc_vec = new long[seq_len];

/***
	printf("\nBeginning The data is: Y_VEC\n");
	for(t = 0; t < 10; t++)
		printf("%lf ,", data.y_vec[t]);
	printf("\nBeginning2 The data is: Y_VEC B\n");
	for(t = 0; t < 10; t++)
		printf("%lf ,", data.y_vecB[t]);
***/

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached 3\n");
	fclose(debug_f);
#endif

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

	in = mxGetPr(prhs[3]);   // Get the  missing data flag
	miss_data = long(in[0]); 

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "Reached 4\n");
	fclose(debug_f);
#endif

	in = mxGetPr(prhs[4]);   // Get the segments starts
	data.num_segments = MAX(mxGetM(prhs[4]), mxGetN(prhs[4]));   // get the number of segments (e.g. chromosomes)

#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	fprintf(debug_f, "NUM SEGMENTS  %ld\n", data.num_segments);
	fclose(debug_f);
#endif
	data.segments_starts = new long[data.num_segments];
	for(t = 0; t < data.num_segments; t++)   
		data.segments_starts[t] = in[t];
/**
	printf("Num Segments : %ld\n", data.num_segments);
	for(t = 0; t < data.num_segments; t++)   
		printf("%ld, ", data.segments_starts[t]);
	printf("\n\n");
**/
#ifdef DEBUG
	debug_f = fopen(debug_p, "a");
	for(t = 0; t < data.num_segments; t++)   
		fprintf(debug_f, "%ld, ", data.segments_starts[t]);
	fclose(debug_f);
#endif

	in = mxGetPr(prhs[5]);   // Get the x-dimension 
	x_dim = long(in[0]);     // Copy the x-dimension

	in = mxGetPr(prhs[6]);   // Get the y-dimension 
	y_dim = long(in[0]);     // Copy the y-dimension
	

	// here insert the place 
	in = mxGetPr(prhs[7]);   // Get the place flag 
	place_flag = long(in[0]);     // Copy the place flag

	// here insert the fold change 
	in = mxGetPr(prhs[8]);				// Get the fold-change flag 
	fold_change_flag = long(in[0]);     // Copy the fold-change flag
	
	in = mxGetPr(prhs[12]);   // Get the special models flag
	hmm.special_models_flag = long(in[0]);


	// Decide if to use x_dim or 36
	long max_x_states = x_dim;
	if(hmm.special_models_flag == SNPS_ONE_SAMPLE)
		max_x_states = 36;



//	printf("hmm.x_dim2 is : %ld  y_dim is %ld max_x_states is %ld seq_len is %ld \n", 
//		hmm.x_dim2, y_dim, max_x_states, seq_len);

//	return;

	if(place_flag)
	{
		hmm.x_dim2 = 2; // just to make sure ..

		// allocate memory 
		if(hmm.special_models_flag == 0)
		{
			for(i = 0; i < max_x_states; i++)		
				for(j = 0; j < y_dim; j++)
				{
					hmm.place_gaussian_mu[i][j] = new double[seq_len];
					hmm.place_gaussian_sigma[i][j] = new double[seq_len];
				}

			// Copy the place means
			in = mxGetPr(prhs[9]);   // Get the place means
			long arr_size = mxGetM(prhs[9]) * mxGetN(prhs[9]);
			for(i = 0; i < max_x_states; i++)
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

			// Copy the place sigmas
			in = mxGetPr(prhs[10]);   // Get the place sigmas
			arr_size = mxGetM(prhs[10]) * mxGetN(prhs[10]);

			for(i = 0; i < max_x_states; i++)	
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
		}

		// Read the place_M
		for(i = 0; i < 2; i++)		
			for(j = 0; j < 2; j++)
			{
				hmm.place_M[i][j] = new double[seq_len];  // conditional probability matrices
				hmm.place_M_cum[i][j] = new double[seq_len]; // cumulative sum
			}
		in = mxGetPr(prhs[11]);   // Get the place M's
		for(i = 0; i < hmm.x_dim2; i++) // here we take the smaller dimension (x_dim2 = 2)	
			for(j = 0; j < hmm.x_dim2; j++)
				for(t = 0; t < seq_len; t++)
					hmm.place_M[i][j][t] = in[t*hmm.x_dim2*hmm.x_dim2+j*hmm.x_dim2+i]; // in[t+seq_len*j+seq_len*hmm.x_dim2*i];
 /***/
	}
	/////////////////
		
//	return;


	// Now determine the upperbounds
	if(hmm.special_models_flag == 0)
	{
		in = mxGetPr(prhs[13]);   // Get the upperbounds
		for(i = 0; i < max_x_states; i++)	
			for(j = 0; j < x_dim; j++)
				hmm.M_upperbounds[i][j] = in[i*x_dim+j];
	}
	in = mxGetPr(prhs[14]);   // Get the upperbounds flag
	hmm.use_bounds = long(in[0]);
	in = mxGetPr(prhs[15]);   // Get the iterations 
	max_iters = long(in[0]); // Copy the iteations
	in = mxGetPr(prhs[16]);   // Get the number of starting points 
	num_starts = long(in[0]); // Copy the number of starting points
	in = mxGetPr(prhs[17]);   // Get the tolerance for stopping iterations 
	tolerance = in[0];		 // Copy the tolerance for stopping iterations
 

	data.seq_len = seq_len; // this does in the data
	data.y_type = CONTINUOUS;
	data.miss_data = miss_data;
	data.dont_miss_prob = 0.5; // not important (only for simulations)



//	printf("X DIM : %ld  Y DIM %ld\n", hmm.x_dim, hmm.y_dim);

	hmm.seq_len = seq_len; 
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
//  	printf("BEFORE EM X DIM : %ld  Y DIM %ld YTYPE: %ld\n", hmm.x_dim, hmm.y_dim, hmm.y_type);

/***
	printf("\nThe data is: Y_VEC\n");
	for(t = 0; t < 10; t++)
		printf("%lf ,", data.y_vec[t]);
	printf("\nThe data is: Y_VEC B\n");
	for(t = 0; t < 10; t++)
		printf("%lf ,", data.y_vecB[t]);
***/
	
//	printf("Calling EM !!! USE BOUNDS : %ld \n", hmm.use_bounds);
//	return; // We don't want to go into the real work yet ..

	// Call c function to do the job for you ..
	if(hmm.special_models_flag == SNPS_ONE_SAMPLE)
	{	
		double *tmp_vec = new double[data.seq_len];
		for(i=0; i<data.seq_len; i++)
			tmp_vec[i] = 0;
		/***
		printf("The place M's:\n");
		for(i=0; i<100; i++)
			printf("%lf %lf %lf %lf | %lf %lf %lf %lf || SUM: %lf\n", hmm.place_M[0][0][i], hmm.place_M[0][1][i],
			hmm.place_M[1][0][i], hmm.place_M[1][1][i],
			hmm.place_M_cum[0][0][i], hmm.place_M_cum[0][1][i],
			hmm.place_M_cum[1][0][i], hmm.place_M_cum[1][1][i], 
			hmm.place_M[0][0][i] + hmm.place_M[0][1][i] + hmm.place_M[1][0][i] + hmm.place_M[1][1][i]);
		/***/

//		printf("DATA LENGTH IS %ld\n", data.seq_len);
		double best_model_index = TrainModelEMSNPs(&hmm, &data, max_iters, num_starts, tolerance, 
			&best_model_score, tmp_vec);

		/***
		printf("AFFFFTER The place M's:\n");
		for(i=0; i<100; i++)
			printf("%lf %lf %lf %lf | %lf %lf %lf %lf || SUM: %lf\n", hmm.place_M[0][0][i], hmm.place_M[0][1][i],
			hmm.place_M[1][0][i], hmm.place_M[1][1][i],
			hmm.place_M_cum[0][0][i], hmm.place_M_cum[0][1][i],
			hmm.place_M_cum[1][0][i], hmm.place_M_cum[1][1][i], 
			hmm.place_M[0][0][i] + hmm.place_M[0][1][i] + hmm.place_M[1][0][i] + hmm.place_M[1][1][i]);
		/***/

/***		
		printf("Used special EM, Best Index is %lf\n", best_model_index); /***
		printf("We got the following scores vec:\n");
		for(i = 0; i < 100; i++)
			printf("%lf, ", tmp_vec[i]); 
		printf("\n"); /***/

/***
		printf("\n The First Scores : \n");
		for(i = 12; i < 30; i++)
			printf("%lf, ", tmp_vec[i]);
		printf("\n The Scale Exps : \n");
		for(i = 30; i < MIN(100,seq_len); i++)
			printf("%lf, ", tmp_vec[i]);
		printf("\n");
/***

		printf("\nThe MU and SIGMA in all iterations\n");
		for(i=0; i < 10; i++)
		{
			for(j=0; j<5; j++)
				printf("%lf ", tmp_vec[i*10+j]);
			printf("  |  ");
			for(j=0; j<5; j++)
				printf("%lf ", tmp_vec[i*10+5+j]);
			printf("\n");
		}
/***/

		delete tmp_vec;
	}
	else
	{
//		printf("Using standard EM DIM is %ld \n", hmm.x_dim);
		TrainModelEM(&hmm, &data, max_iters, num_starts, tolerance, &best_model_score);
	}

	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/***	printf("HMM MODEL : \n");
	printf("X DIM : %ld  Y DIM %ld YTYPE: %ld\n", hmm.x_dim, hmm.y_dim, hmm.y_type);

	hmm.x_dim = 3; hmm.y_dim = 1; // force dimensions back for debugging	
	if(hmm.y_type == DISCRETE)
	{
		printf("Y DIM : %ld\n", hmm.y_dim);
		printf("X is Discrete, Y is Discrete\n");
	}
	else
		printf("X is Discrete, Y is Mixture of %ld Gaussians\n", hmm.y_dim);
	printf("Initial X distribution : \n");
	for(i = 0; i < hmm.x_dim; i++)
		printf("%lf, ", hmm.PI[i]);
	printf("\n");

	printf("X Transition Matrix :\n");
	for(i = 0; i < hmm.x_dim; i++)
	{
		for(j = 0; j < hmm.x_dim; j++)
			printf("%lf, ", hmm.M[i][j]);
		printf("\n");
	}

	printf("X CUMCUM Transition Matrix :\n");
	for(i = 0; i < hmm.x_dim; i++)
	{
		for(j = 0; j < hmm.x_dim; j++)
			printf("%lf, ", hmm.M_cum[i][j]);
		printf("\n");
	}
	if(hmm.y_type == DISCRETE)
	{
		printf("X->Y Omission Matrix :\n");
		for(i = 0; i < hmm.x_dim; i++)
			PrintDoubleVec(hmm.N[i], hmm.y_dim);
		printf("X->Y CUMCUM Omission Matrix :\n");
		for(i = 0; i < hmm.x_dim; i++)
			PrintDoubleVec(hmm.N_cum[i], hmm.y_dim);
	}
	else
	{
		printf("X->Y Mixture Matrix :\n");
		for(i = 0; i < 2*hmm.x_dim-1; i++)
		{
			for(j = 0; j < hmm.y_dim; j++)
				printf("%lf, ",hmm.N[i][j]);
	//		printf("\n");
		}
		printf("\nY Mean Matrix :\n");
		for(i = 0; i < 2*hmm.x_dim-1; i++)
		{
			for(j = 0; j < hmm.y_dim; j++)
				printf("%lf, ",hmm.MU[i][j]);
//			printf("\n");
		}
		printf("\nY Std. Matrix :\n");
		for(i = 0; i < 2*hmm.x_dim-1; i++)
		{
			for(j = 0; j < hmm.y_dim; j++)
				printf("%lf, ",hmm.SIGMA[i][j]);
//			printf("\n");
		}

	}
****/   // ended printing
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	long eq_dim = hmm.x_dim +  hmm.special_models_flag * (hmm.x_dim-1);

	// Now copy output from the HMM structure  : 
	plhs[0] = mxCreateDoubleMatrix(1, x_dim, mxREAL);  // The initial vector PI
	out = mxGetPr(plhs[0]); // output PI
	for(i=0; i<hmm.x_dim; i++) 
		out[i] = hmm.PI[i];
	plhs[1] = mxCreateDoubleMatrix(x_dim, x_dim, mxREAL);  // The transition matrix M
	out = mxGetPr(plhs[1]); // output M
	for(i=0; i<x_dim; i++) 
		for(j=0; j<x_dim; j++)
			out[j*x_dim+i] = hmm.M[i][j];
	plhs[2] = mxCreateDoubleMatrix(eq_dim, y_dim, mxREAL);  // The mixture matrix N
	out = mxGetPr(plhs[2]); // output N
	for(i = 0; i < eq_dim; i++) 
		for(j = 0; j < y_dim; j++)
			out[j*eq_dim+i] = hmm.N[i][j];
	plhs[3] = mxCreateDoubleMatrix(eq_dim, y_dim, mxREAL);  // The mean matrix MU
	out = mxGetPr(plhs[3]); // output MU
	for(i = 0; i < eq_dim; i++) 
		for(j = 0; j < y_dim; j++)
			out[j*eq_dim+i] = hmm.MU[i][j];
	plhs[4] = mxCreateDoubleMatrix(eq_dim, y_dim, mxREAL);  // The mixture matrix SIGMA
	out = mxGetPr(plhs[4]); // output SIGMA
	for(i = 0; i < eq_dim; i++) 
		for(j = 0; j < y_dim; j++)
			out[j*eq_dim+i] = hmm.SIGMA[i][j];
	plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);  // The log-score
	out = mxGetPr(plhs[5]); // output log-score
	out[0] = best_model_score;


	// free memory
	delete data.y_vec;
	delete data.y_vecB;
	delete data.loc_vec; 
	delete data.loc_diff_vec;
	delete data.segments_starts;

	if(place_flag)
	{	

		// delete memory 
		if(hmm.special_models_flag == 0)
			for(i = 0; i < hmm.x_dim; i++)		
				for(j = 0; j < hmm.y_dim; j++)
				{
					delete hmm.place_gaussian_mu[i][j];
					delete hmm.place_gaussian_sigma[i][j];
				}
		for(i = 0; i < hmm.x_dim2; i++)		
			for(j = 0; j < hmm.x_dim2; j++)
			{
				delete hmm.place_M[i][j]; 
				delete hmm.place_M_cum[i][j]; 
			}

	}


}
//#endif // MEX_TO_MATLAB




