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
  
	int nRows, nCols, i, j, K, t, max_iters, num_starts;

	long seq_len, L; 
	long x_dim, y_dim, place_flag;
	double tolerance;
  
	double *in, *out;

	long *y_vec_int; // dummy !!! 
	double *y_vec; 
	double *loc_vec;  // location on chromosome vector - Currently not used !!!! 
	
	hmm_model hmm;

	double *b_cond_tabs[MAX_X_VALS]; // helpful conditional tables
	double *alpha[MAX_X_VALS];
	double *beta[MAX_X_VALS];
	double *gamma[MAX_X_VALS];
	double *gammagamma[MAX_X_VALS][MAX_Y_VALS];
	double *y_mew_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS];


	long *opt_x_vec_out;  // optimal viterbi path

	double *scale, *scale_cum, *scale_exp;
	double cur_score;

	 /* Check for proper number of arguments. Should be Seven 7 */
	if(nrhs != 10) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. Expression vector
	2. Chromosomal location vector
	3. PI initial vector
	4. M transition matrix
	5. N mixture matrix 
	6. MEW mean matrix
	7. SIGMA standard error matrix
	8. place flag
	9. place means
	10. place sigmas

 ---   Note : Dimensions are already IN the model !!! ---
	*/


	in = mxGetPr(prhs[0]);   // Get the expression vector
	seq_len = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the expression vector length

	y_vec = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy expression vector
		y_vec[t] = in[t];  
	


	in = mxGetPr(prhs[1]);   // Get the location vector
	L  = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));    // Get the number of sequences (genes)
	if(L != seq_len)
		mexErrMsgTxt("Usage: Exp. and Loc. vector length must be the same");

	loc_vec = new double[seq_len];
	for(t = 0; t < seq_len; t++)   // copy location vector
		loc_vec[t] = in[t];  

	

	in = mxGetPr(prhs[2]);   // Get the initial vector
	x_dim = MAX(mxGetM(prhs[2]), mxGetN(prhs[2]));     // Get the x dim

  	for(i = 0; i < x_dim; i++)
		hmm.PI[i] = in[i];   // Copy the initial vector


	in = mxGetPr(prhs[3]);   // Get the transition matrix
	x_dim = mxGetM(prhs[3]);     // Get the x dim
	L = mxGetN(prhs[3]);
	if(L != x_dim)
		mexErrMsgTxt("Usage: transition matrix must be square ");

  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < x_dim; j++)
			hmm.M[i][j] = in[j*x_dim+i];   // Copy the transition matrix


	in = mxGetPr(prhs[4]);   // Get the mixture matrix
	L = mxGetM(prhs[4]);
	y_dim = mxGetN(prhs[4]);     // Get the y dim
	if(L != x_dim)
		mexErrMsgTxt("Usage: emission and transition matrices must have same number of rows ");
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.N[i][j] = in[j*x_dim+i];   // Copy the mixture matrix

	in = mxGetPr(prhs[5]);   // Get the mean matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.MEW[i][j] = in[j*x_dim+i];   // Copy the mean matrix

	in = mxGetPr(prhs[6]);   // Get the std matrix
  	for(i = 0; i < x_dim; i++)
		for(j = 0; j < y_dim; j++)
			hmm.SIGMA[i][j] = in[j*x_dim+i];   // Copy the std matrix

	in = mxGetPr(prhs[7]);   // Get the place flag
	place_flag = mxGetM(prhs[7]);


	printf("BEFORE PLACE FLAG \n");

	if(place_flag)
	{
		
		printf("INSIDE PLACE FLAG \n");
		// allocate memory 
		for(i = 0; i < x_dim; i++)		
			for(j = 0; j < y_dim; j++)
			{
				printf("alloc (%ld %ld)\n", i, j);
				hmm.place_gaussian_mew[i][j] = new double[seq_len];
				hmm.place_gaussian_sigma[i][j] = new double[seq_len];
			}


		// Now copy
		in = mxGetPr(prhs[5]);   // Get the place means
		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
				for(t = 0; t < seq_len; t++)
					hmm.place_gaussian_mew[i][j][t] = in[0]; ///in[t+seq_len*j+seq_len*hmm.y_dim*i]; // Copy the place means


		// Now copy
		in = mxGetPr(prhs[6]);   // Get the place sigmas
		for(i = 0; i < x_dim; i++)	
			for(j = 0; j < y_dim; j++)
			{
				printf("(I J ) %ld %ld \n", i, j);
				for(t = 0; t < seq_len; t++)
					hmm.place_gaussian_sigma[i][j][t] = in[0]; ///in[t+seq_len*j+seq_len*hmm.y_dim*i]; // Copy the place sigmas
			}


	}

	printf("AFTER PLACE FLAG \n");



	hmm.seq_len = seq_len;
	hmm.x_dim = x_dim;
	hmm.y_dim = y_dim;
	hmm.y_type = CONTINUOUS;
	hmm.cum_flag = 1;
	hmm.log_flag = 1;
	hmm.place_flag = place_flag; // currently no placing 



	ComputeModelAuxillaryParameters(&hmm); // Compute the cums ...

	// allocate memory
	for(i = 0; i < x_dim; i++)
	{
		b_cond_tabs[i] = new double[seq_len];
		alpha[i] = new double[seq_len];
		beta[i] = new double[seq_len];
		gamma[i] = new double[seq_len];

		for(j = 0; j < y_dim; j++)
		{
			gammagamma[i][j] = new double[seq_len];
			y_mew_square_exp_vecs[i][j] = new double[seq_len];
		}
	}
	opt_x_vec_out = new long[seq_len];
	scale = new double[seq_len];
	scale_cum = new double[seq_len];
	scale_exp = new double[seq_len];

	// Call many c functions to do the job for you ..

	Compute_b_cond_tabs(&hmm, y_vec, y_vec_int, b_cond_tabs, y_mew_square_exp_vecs);  // OK

	

	// Now find the probs
	forward(&hmm, b_cond_tabs, alpha, scale, scale_cum, scale_exp, &cur_score);  // PROBLEM FIRST ALPHA !!! 
	backward(&hmm, scale_exp, b_cond_tabs, beta);		

	ComputeMarginalGamma(&hmm, y_vec, y_mew_square_exp_vecs, alpha, beta, scale_exp, b_cond_tabs,
						 gamma, gammagamma);

	/** Print some ...***/
//	for(i = 0;i < x_dim; i++)
//		for(t=0; t<10; t++) 
//			printf("Gamma[%ld][%ld] %lf alph  %lf bet %lf b_con %lf scale %lf\n", 
//			i, t, gamma[i][t], alpha[i][t], beta[i][t], b_cond_tabs[i][t], scale[t]);


//	printf("Inside Score : %lf\n", cur_score);

//	printf("Start Viterbi\n");
	Viterbi( &hmm, b_cond_tabs, opt_x_vec_out);

//	printf("Start COPY OUTPUT\n");
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
		delete b_cond_tabs[i];
		delete alpha[i];
		delete beta[i];
		delete gamma[i];
		for(j = 0; j < y_dim; j++)
		{
			delete gammagamma[i][j];
			delete y_mew_square_exp_vecs[i][j];
		}
	}
	delete opt_x_vec_out;
	delete y_vec;
	delete loc_vec; 
	delete scale;
	delete scale_cum; 
	delete scale_exp;
		

	if(place_flag)
	{	

		// delete memory 
		for(i = 0; i < hmm.x_dim; i++)		
			for(j = 0; j < hmm.y_dim; j++)
			{
				delete hmm.place_gaussian_mew[i][j];
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
			PrintDoubleVec(hmm->MEW[i], hmm->y_dim);
		printf("Y Std. Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec(hmm->SIGMA[i], hmm->y_dim);
	}



	return 0;

}











#ifdef EVERYTHING_INSIDE




#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "general.h"
#include "hmm_chrom_funcs.h"












/**************************************************************************/

/******************************
 * qsort.c                    *
 * quicksort algorithm in C   *
 * Daniel J. Schultz          *
 ******************************/

// inner split routine
long split(double *vals, long first, long last, long *indexes)
{
    double    x, temp;
    long      u, sp, temp_int;

    x  = vals[first];
    sp = first;
    for (u = first + 1; u <= last; ++u)
       if( vals[u] < x)
	   {
           sp++;
           SWAP(vals[u], vals[sp], temp);
		   SWAP(indexes[u], indexes[sp], temp_int);
       }

    SWAP(vals[first], vals[sp], temp);
	SWAP(indexes[first], indexes[sp], temp_int);
    
	return(sp);
}

// recursive quicksort routine
long quicksort(double *vals, long first, long last, long *indexes)
{
    long     splitpt;

    if (first < last) 
	{
       splitpt = split(vals, first, last, indexes);
       quicksort(vals, first, splitpt - 1, indexes);
       quicksort(vals, splitpt +1, last, indexes);
    }

	return 0;
}

// outer  call for quicksort
long DoQuicksort(double *vals, long len, long *indexes)
{
	long i;

	for(i = 0; i < len; i++)
		indexes[i] = i;  // init the indexes 
		
	quicksort(vals, 0, len-1, indexes);

	return 0; 
}










long PrintVec(long *vec, long len)
{
	long i;
	for(i = 0; i < len; i++)
	{
		printf("%ld ", vec[i]);
		if((i%PRINT_IN_ROW) == PRINT_IN_ROW-1)
			printf("\n");
	}

	printf("\n");

	return 0; 
}


long PrintDoubleVec(double *vec, long len)
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


double randr()
{
	return double(rand())/double((RAND_MAX+1.0));
}


// Generates a gaussian random variable with mean mew and standard deviation sigma 
double Gaussian(double mew,double sigma)

{
   double t = 0.0;
   double x,v1,v2,r;
   
   if (t == 0) 
   {
		do 
		{
			v1 = 2.0 * randr() - 1.0;
			v2 = 2.0 * randr() - 1.0;
			r = v1 * v1 + v2 * v2;
		}  while (r >= 1.0);

		r = sqrt((-2.0*log(r))/r);
		t = v2*r;
		return(mew+v1*r*sigma);
   }
   else 
   {
     x = t;
     t = 0.0;
     return(mew+x*sigma);
   }
 }



// Perform the Viterbi algorithm. Find the most probable path for the X's given a sequence of Y's
// Implementation is according to Rabiner's tutorial
// Everything is done with logs to avoid underflows
long Viterbi( hmm_model *hmm, double *b_cond_tabs[MAX_X_VALS], long seq_len, 
			  long *opt_x_vec_out)
{
//	printf("SSSSSStart Viterbi\n");

	long i, j, t;
	
	double *lambda[MAX_X_VALS];
	long *max_indexes[MAX_X_VALS];

	double cur_log_prob, max_log_prob, cur_log_prob_lambda, max_log_prob_lambda;
	long max_index;

//	printf("Start Viterbi\n");


	// First allocate memory 
	for(i = 0; i < hmm->x_dim; i++)
	{
		lambda[i] = new double[seq_len];
		max_indexes[i] = new long[seq_len];
	}
		


	// Initilize
	for(j = 0; j < hmm->x_dim; j++)
	{
		lambda[j][0] = log(hmm->PI[j]) + log(b_cond_tabs[j][0]);		
		max_indexes[j][0] = -1;   // maybe 0/-1 ???
	}

		
//	printf("Done Init Viterbi\n");

	// Perform recursive procedure
	for(t = 1; t < seq_len; t++)
	{
		// Find the maximal lambda

		for(j = 0; j < hmm->x_dim; j++)			// state at time t
		{
			max_index = -1; max_log_prob = -9999999999999; max_log_prob_lambda = -999999999999;
			for(i = 0; i < hmm->x_dim; i++)		// state at time t-1
			{
				cur_log_prob = lambda[i][t-1] + log(hmm->M[i][j]);
				cur_log_prob_lambda = cur_log_prob +  log(b_cond_tabs[j][t]); 
				if(cur_log_prob > max_log_prob)  // find the max
				{
					max_log_prob = cur_log_prob;
					max_index = i;
				}
				if(cur_log_prob_lambda > max_log_prob_lambda)  // find the max
					max_log_prob_lambda = cur_log_prob_lambda;

			}
			lambda[j][t] = max_log_prob_lambda;  // take the maximal prob.
			max_indexes[j][t] = max_index;   // take the maximal index
		}
	}

//	printf("Done all letters Viterbi\n");


//	for(i = 0; i < hmm->x_dim; i++)
//	{
//		print_double_vec(lambda[i], seq_len);
//		print_vec(max_indexes[i], seq_len);
//	}


	// Termination 
	double p_max = -99999999999999.9;
	long p_max_index = -1;
	for(i = 0; i < hmm->x_dim; i++)
		if(p_max < lambda[i][t-1])
		{
	//		printf("FOUND PMAX !\n");
			p_max = lambda[i][t-1];
			p_max_index = i;
		}


//	printf("Done pmax Viterbi\n");

	// Find optimal path using backtracking
	opt_x_vec_out[seq_len-1] = p_max_index;
	
	for(t = seq_len-2; t >= 0; t--)
	{
//		printf("ttt %ld opttt %ld \n", t, opt_x_vec_out[t+1]);
		opt_x_vec_out[t] = max_indexes[opt_x_vec_out[t+1]][t+1];
	}
//	printf("Done backtrack Viterbi\n");

	// Last free memory 
	for(i = 0; i < hmm->x_dim; i++)
	{
		delete lambda[i];
		delete max_indexes[i];
	}


	return 0; 

}



// Perform the Forward algorithm. Find the alpha coefficients and the prob. of Y (summing over all possible X's)
// Implementation is according to Rabiner's tutorial
// need to do 'scaling' of the alphas to avoid underflows ..
long forward( hmm_model *hmm, long seq_len, double *b_cond_tabs[MAX_X_VALS], 
			  double *alpha[MAX_X_VALS], double *scale, double *scale_cum,  double *scale_exp, double *y_vec_log_prob_out)
{
	long i, j, t;
	

	// Initilize
	for(j = 0; j < hmm->x_dim; j++)
		alpha[j][0] = hmm->PI[j] * b_cond_tabs[j][0];    // use the already computed b tables 
	scale[0] = 0; scale_exp[0] = 1;


	// Perform recursive procedure	
	for(t = 1; t < seq_len; t++)
	{
		// Compute the next alphas. Note : We need scaling here !!!
		if((t%DO_SCALE)==DO_SCALE-1)
		{
			scale_exp[t]=0;
			for(j = 0; j < hmm->x_dim; j++)			// state at time t
			{
				alpha[j][t] = alpha[0][t-1] * hmm->M[0][j]; // take the 1st one ..	
				for(i = 1; i < hmm->x_dim; i++)		// state at time t-1
					alpha[j][t] += alpha[i][t-1] * hmm->M[i][j];
				alpha[j][t] *= b_cond_tabs[j][t];
				scale_exp[t] += alpha[j][t];
			}
		
			// Now do scaling	
			for(j = 0; j < hmm->x_dim; j++)
				alpha[j][t] /= scale_exp[t];
			scale[t] = log(scale_exp[t]); // transfer to logs ..
		}
		else  // No Scailing !!
		{
			scale[t] = 0; scale_exp[t] = 1;
			for(j = 0; j < hmm->x_dim; j++)			// state at time t
			{
				alpha[j][t] = alpha[0][t-1] * hmm->M[0][j]; // take the 1st one ..	
	
				for(i = 1; i < hmm->x_dim; i++)		// state at time t-1
					alpha[j][t] += alpha[i][t-1] * hmm->M[i][j];
				alpha[j][t] *= b_cond_tabs[j][t];
			}
		}
	}


	// Termination 
	*y_vec_log_prob_out = 0;


	// scaled alpha's : 
	scale_cum[0] = scale[0];
	for(t = 1; t < seq_len; t++)
		scale_cum[t] = scale_cum[t-1] + scale[t];

	*y_vec_log_prob_out = scale_cum[seq_len-1];

/*
	if(scale_cum[seq_len-1] > 0)  // Positive log-likelihood !!! Impossible !!!! 
	{
		printf("SCALE : \n");
		PrintDoubleVec(scale, seq_len);
		printf("SCALE CUM: \n");
		PrintDoubleVec(scale_cum, seq_len);
		printf("SCALE EXP : \n");
		PrintDoubleVec(scale_exp, seq_len);
	}
*/

	return 0; 

}




// Perform the Backward algorithm. Find the beta coefficients and the prob. of Y (summing over all possible X's)
// Implementation is according to Rabiner's tutorial
// Note : We assume here that the 'scaling' coefficients are already determined by the forward algorithm
long backward( hmm_model *hmm, long seq_len, double *scale_exp, double *b_cond_tabs[MAX_X_VALS],
			  double *beta[MAX_X_VALS])
{
	long i, j, t;
	

//	printf("Start back\n");

	// Initilize
	for(i = 0; i < hmm->x_dim; i++)
		beta[i][seq_len-1] = 1 / scale_exp[seq_len-1]; // Here we start with one scaled.		
	
	for(t = seq_len-2; t >= 0; t--)      // Compute the next betas. Note : We need scaling here !!!
		for(i = 0; i < hmm->x_dim; i++)			// state at time t
		{
			beta[i][t] = beta[0][t+1] * hmm->M[i][0] *  b_cond_tabs[0][t+1]; // take the 1st one ..
			for(j = 1; j < hmm->x_dim; j++)		// state at time t+1
				beta[i][t] += beta[j][t+1] * hmm->M[i][j] * b_cond_tabs[j][t+1]; 
			beta[i][t] /= scale_exp[t];  // already apply scaling !!!
		}

//	printf("Finish back\n");


	// Termination - This is only for testing !!!
//	*y_vec_log_prob_out = 0;
//
//	// scaled beta's : we compute again the prob. just for fun ..
//	for(t = 0; t < seq_len; t++)
//		*y_vec_log_prob_out += scale[t];



	return 0; 

}


// Here we find thge marginal prob. of each X variable, This is an alternative to the Viterbi algorithm,
// Note that if interpreted incorrectly it can give 'illegal path' !
// We assume here we already done the forward-backward procedure !!
long ComputeMarginalGamma(hmm_model *hmm,  double *y_vec, long seq_len, double *y_mew_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS],
						 double *alpha[MAX_X_VALS], double *beta[MAX_X_VALS], double *scale_exp, double *b_cond_tabs[MAX_X_VALS],
						 double *gamma[MAX_X_VALS], double *gammagamma[MAX_X_VALS][MAX_Y_VALS])
{
	long i, j, t;

	for(i = 0; i < hmm->x_dim; i++)
		for(t = 0; t < seq_len; t++)
			gamma[i][t] = alpha[i][t] * beta[i][t] * scale_exp[t]; // here scale is already in exp !!! // y_vec_log_prob


	// here we compute gamma gamma !!  we use normal density here !!
	if(hmm->y_type == CONTINUOUS)
	{
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				for(t = 0; t < seq_len; t++)
				{
			//    Removed for performance. Put back if there're problems !!!! 
			//		if(gamma[i][t] == 0)
			//			gammagamma[i][j][t] = 0;
			//		else



					{
				/**	gammagamma[i][j][t] = gamma[i][t] * hmm->N[i][j] *  ONE_OVER_2SQRT_PI * 
					exp (- (y_vec[t]-hmm->MEW[i][j])*(y_vec[t]-hmm->MEW[i][j])/(2*hmm->SIGMA[i][j]*hmm->SIGMA[i][j])) / 
					( hmm->SIGMA[i][j] * b_cond_tabs[i][t]);

					
				/**	gammagamma[i][j][t] = gamma[i][t] * hmm->N[i][j] *  ONE_OVER_2SQRT_PI * 
					exp (- (y_vec[t]-hmm->MEW[i][j])*(y_vec[t]-hmm->MEW[i][j])*0.5*hmm->SIGMA_inv[i][j]*hmm->SIGMA_inv[i][j]) / 
					( hmm->SIGMA[i][j] * b_cond_tabs[i][t]);
				**/	

////////////////////////////////////////
					gammagamma[i][j][t] = gamma[i][t] * ///// ALREADY hmm->N[i][j] *  ONE_OVER_2SQRT_PI * 
					y_mew_square_exp_vecs[i][j][t] / b_cond_tabs[i][t]; 
///// ALREADY					( hmm->SIGMA[i][j] * b_cond_tabs[i][t]);
////////////////////////////////////////					
					
					
  
					}
				}

	/**	
//		double sumsum;
//			// Sanity check !!!!
//		for(i = 1 /*0*; i < 2/*hmm->x_dim*; i++)
//			for(t = 15 /*0*; t < 20/*seq_len*; t++)
//			{
//				printf("t %ld GAMMA_DIRECT %lf ", t, gamma[i][t]);
//				sumsum = 0;
//
//				for(j = 0; j < hmm->y_dim; j++)
//					sumsum += gammagamma[i][j][t];
//				printf("GAMMA_SUM %lf  DIFF %lf \n", sumsum, sumsum-gamma[i][t]);
//			}
//
//	 /**/
	}


	
	return 0;

}


// Compute the auxillary psi parameters used to re-estimate the model parameters later
long ComputePsi(hmm_model *hmm, long seq_len, double *alpha[MAX_X_VALS], double *beta[MAX_X_VALS], double *scale, double *b_cond_tabs[MAX_X_VALS],
				 double *psi[MAX_X_VALS][MAX_X_VALS])
{
	long i, j, t;

	double temp;

	// a simple loop ...
	for(i = 0; i < hmm->x_dim; i++)
		for(j = 0; j < hmm->x_dim; j++)
			for(t = 0; t < seq_len-1; t++)		
				psi[i][j][t] = hmm->M[i][j] * alpha[i][t] * beta[j][t+1] * b_cond_tabs[j][t+1];

	return 0; 

}


// This is only done for the continunuous case : Here we replace the tabs of N by computing
// Mixture of Gaussians coefficients !! 
// Here b_cond_tabs[i][t] = Pr(Y = y_vec[t] | X = i = \sum_{j}  N[i][j] * Gaussian_{i,j} (y_vec[t]) )
// Note : If we want tommorow a different continuous distribution, we chagne only (?) here !!! 
// We compute here also for discrete distributions to gain (?) performance !!! 
long Compute_b_cond_tabs(hmm_model *hmm, double *y_vec, long *y_vec_int, long seq_len, 
					double *b_cond_tabs[MAX_X_VALS], double *y_mew_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS])
{
	long i, j, t;

	if(hmm->y_type == DISCRETE)
	{
		for(i = 0; i < hmm->x_dim; i++)
			for(t = 0; t < seq_len; t++)
				b_cond_tabs[i][t] = hmm->N[i][y_vec_int[t]];
	}
	else
	{
/**/

		if(hmm->place_flag == 0)  // standard model. No problem ...	
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					for(t = 0; t < seq_len; t++)
						y_mew_square_exp_vecs[i][j][t] =  ONE_OVER_2SQRT_PI * hmm->N[i][j] * hmm->SIGMA_inv[i][j] *
						exp( -(y_vec[t] - hmm->MEW[i][j])*(y_vec[t] - hmm->MEW[i][j]) * 0.5 *hmm->SIGMA_inv[i][j]*hmm->SIGMA_inv[i][j] ); 
		else   // here different gaussian for every place 
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					for(t = 0; t < seq_len; t++)
						y_mew_square_exp_vecs[i][j][t] =  ONE_OVER_2SQRT_PI * hmm->N[i][j] * 
						exp( -(y_vec[t] - hmm->place_gaussian_mew[i][j][t])*(y_vec[t] - hmm->place_gaussian_mew[i][j][t]) / 
						(2.0*hmm->place_gaussian_sigma[i][j][t]*hmm->place_gaussian_sigma[i][j][t]) ) / hmm->place_gaussian_sigma[i][j][t]; 

/**/



		for(i = 0; i < hmm->x_dim; i++)
		{
			for(t = 0; t < seq_len; t++)
				b_cond_tabs[i][t] = 0;
			for(j = 0; j < hmm->y_dim; j++)
				for(t = 0; t < seq_len; t++)
				{

				

					{

					

/*					b_cond_tabs[i][t] += 
					hmm->N[i][j] * ONE_OVER_2SQRT_PI * 
					exp( -(y_vec[t] - hmm->MEW[i][j])*(y_vec[t] - hmm->MEW[i][j]) / (2.0*hmm->SIGMA[i][j]*hmm->SIGMA[i][j]) ) / 
					hmm->SIGMA[i][j];

*/
						b_cond_tabs[i][t] += 
					//	hmm->N[i][j] * 
						y_mew_square_exp_vecs[i][j][t];/// *
					//	exp( -y_mew_square_exp_vecs[i][j][t] * 
					//	hmm->SIGMA_inv[i][j]*hmm->SIGMA_inv[i][j]*0.5 ) * 
					//	hmm->SIGMA_inv[i][j];
					

					}						
	
				//////	b_cond_tabs[i][t] *= ONE_OVER_2SQRT_PI;
		//					exp( hmm->N_log[i][j] + ONE_OVER_2SQRT_PI_LOG + 
		//					 (y_vec[t] - hmm->MEW[i][j])*(y_vec[t] - hmm->MEW[i][j]) / (2.0*hmm->SIGMA[i][j]*hmm->SIGMA[i][j])  -
		//					hmm->SIGMA_log[i][j] );     // add the mixture coefficient
				}

	////		for(t = 0; t < seq_len; t++)
	////			b_cond_tabs[i][t] *= ONE_OVER_2SQRT_PI;
	//// ALREADY IN THE EXPONENT Y_EXP !!! 

		}

	}


	// Now make sure that they are not too small !!!!!
	for(i = 0; i < hmm->x_dim; i++)
		for(t = 0; t < seq_len; t++)
			b_cond_tabs[i][t] = MAX(b_cond_tabs[i][t], EPSILON);
		


	return 0; 

}


// Update the current model parameters to improve likelihood !
long UpdateModelParams(hmm_model *hmm, double *y_vec, long *y_vec_int, long seq_len, double *y_mew_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS],
						 double *gamma[MAX_X_VALS], double *gammagamma[MAX_X_VALS][MAX_Y_VALS],
						 double *psi[MAX_X_VALS][MAX_X_VALS])
{
	long i, j, t;


	double gamma_sum[MAX_X_VALS];
	double psi_sum[MAX_X_VALS][MAX_X_VALS];
	double gamma_sum_for_given_output[MAX_X_VALS][MAX_Y_VALS];



	// Compute the sums for the gammas and psis. Note : We go only until T-1 !!
	for(i = 0; i < hmm->x_dim; i++)
	{
		// Init to zero 
		gamma_sum[i] = 0;
		for(j = 0; j < hmm->x_dim; j++)
			psi_sum[i][j] = 0; 

		// loop to sum
		for(t = 0; t < seq_len-1; t++)
		{
			gamma_sum[i] += gamma[i][t];
			for(j = 0; j < hmm->x_dim; j++)
				psi_sum[i][j] += psi[i][j][t];
		}

		if(hmm->y_type == DISCRETE)
		{
			for(j = 0; j < hmm->y_dim; j++)
				gamma_sum_for_given_output[i][j] = 0; 
			for(t = 0; t < seq_len-1; t++)
				gamma_sum_for_given_output[i][y_vec_int[t]] += gamma[i][t];
		}

		/*
		for(j = 0; j < hmm->x_dim; j++)
		{
			psi_sum[i][j] = 0; 
			for(t = 0; t < seq_len-1; t++)
				psi_sum[i][j] += psi[i][j][t];			
		}
		*/
	}

	// First update the initial distribution PI
	for(i = 0; i < hmm->x_dim; i++)
		hmm->PI[i] = gamma[i][0];



	if(hmm->update_M_flag) 
	{
		// Now update the transition probability matrix M
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				hmm->M[i][j] = psi_sum[i][j] / gamma_sum[i];


		// avoid too little probs !!!!
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				hmm->M[i][j] = MAX(EPSILON, hmm->M[i][j]);


	}


	// Finally update the emission matrix N
	for(i = 0; i < hmm->x_dim; i++)
		gamma_sum[i] += gamma[i][seq_len-1];

	if(hmm->y_type == DISCRETE)
	{

		if(hmm->update_N_flag)
		{
			// add the last 
			for(i = 0; i < hmm->x_dim; i++)
				gamma_sum_for_given_output[i][y_vec_int[seq_len-1]] += gamma[i][seq_len-1];

			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->N[i][j] = gamma_sum_for_given_output[i][j] / gamma_sum[i];
		}
	}
	else  // CONTINUOUS
	{
		// Here the continuous (MoG) case : 

		double gammagamma_sum[MAX_X_VALS][MAX_Y_VALS];
		double gammagamma_obs_sum[MAX_X_VALS][MAX_Y_VALS];
		double gammagamma_obs_var_sum[MAX_X_VALS][MAX_Y_VALS];


		// compute gammagamma sum
		if(hmm->place_flag == 0) // gaussians independent of place
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
				{
					gammagamma_sum[i][j] = 0; gammagamma_obs_sum[i][j] = 0; gammagamma_obs_var_sum[i][j] = 0;
					for(t = 0; t < seq_len-1; t++)
					{
						gammagamma_sum[i][j] += gammagamma[i][j][t];
						gammagamma_obs_sum[i][j] += gammagamma[i][j][t] * y_vec[t];
						gammagamma_obs_var_sum[i][j] += gammagamma[i][j][t] * 
							(y_vec[t] - hmm->MEW[i][j]) * (y_vec[t] - hmm->MEW[i][j]);
					}
				}
		else   // here everything depends on place
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
				{
					gammagamma_sum[i][j] = 0; gammagamma_obs_sum[i][j] = 0; gammagamma_obs_var_sum[i][j] = 0;
					for(t = 0; t < seq_len-1; t++)
					{
						gammagamma_sum[i][j] += gammagamma[i][j][t];
						gammagamma_obs_sum[i][j] += gammagamma[i][j][t] * y_vec[t];
						gammagamma_obs_var_sum[i][j] += gammagamma[i][j][t] * 
							(y_vec[t] -  hmm->place_gaussian_mew[i][j][t]) * (y_vec[t] -  hmm->place_gaussian_mew[i][j][t]); 						
					}
				}



		if(hmm->update_N_flag)
		{

			// Update the N's (mixture coefficients)
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->N[i][j] = gammagamma_sum[i][j] / gamma_sum[i];

			// Normalize by sum 
			double sumsum;
			for(i = 0; i < hmm->x_dim; i++)
			{
				sumsum = 0; 
				for(j = 0; j < hmm->y_dim; j++)
					sumsum += hmm->N[i][j];
				for(j = 0; j < hmm->y_dim; j++)
					hmm->N[i][j] /= sumsum;
			}
		}

		
		if(hmm->update_MEW_flag)
			// Update MEW's (the means)
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->MEW[i][j] = gammagamma_obs_sum[i][j]/gammagamma_sum[i][j];

		
		if(hmm->update_SIGMA_flag)
		{
			// Update SIGMA's (the standard deviations)
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->SIGMA[i][j] = sqrt(gammagamma_obs_var_sum[i][j]/gammagamma_sum[i][j]);


			// avoid too little sigma !!!!
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->SIGMA[i][j] = MAX(EPSILON, hmm->SIGMA[i][j]);
		}


	}  // end if continuous


	// avoid too little probs !!!!
	if(hmm->update_N_flag)
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				hmm->N[i][j] = MAX(EPSILON, hmm->N[i][j]);


	ComputeModelAuxillaryParameters(hmm);  // Must be last !!!

	return 0; 
			
}	



// Get a smart starting point for the hmm parameters ..
// Not implemented yet ..
long SmartStartParams(hmm_model *hmm, double *y_vec, long seq_len)
{
	return 0;
}



// perform the EM algorithm to establish the model parameters from data
long TrainModelEM(hmm_model *hmm, double *y_vec, long seq_len, long max_iters, long num_starts, double tolerance)
{
	
	printf("START EM !\n");
	
	long i, j, t;


	// First Initilize model, use the same function so we have here a 'self-reference paradox' ..
	long x_dim = hmm->x_dim;
	long y_dim = hmm->y_dim;
	long y_type = hmm->y_type;
	long place_flag = hmm->place_flag;

	hmm->cum_flag = 0;
	hmm->log_flag = 1;



	// Now start the iterative EM ;-)
	double *b_cond_tabs[MAX_X_VALS];
	double *alpha[MAX_X_VALS];
	double *beta[MAX_X_VALS];
	double *gamma[MAX_X_VALS];
	double *gammagamma[MAX_X_VALS][MAX_Y_VALS];
	double *psi[MAX_X_VALS][MAX_X_VALS];	// The estimated transition probs.
	double *scale = new double[seq_len];   // scaling to avoid underflows
	double *scale_cum = new double[seq_len];   // cumulative sum of the scale
	double *scale_exp = new double[seq_len];   // exp of the scale (ratio)
	double *y_mew_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS];    // (y-mew)^2



	long *y_vec_int = new long[seq_len];
	for(t = 0; t < seq_len; t++)
		y_vec_int[t] = long(y_vec[t]);

	// First allocate memory 
	for(i = 0; i < hmm->x_dim; i++)
	{
		b_cond_tabs[i] = new double[seq_len];
		alpha[i] = new double[seq_len];
		beta[i] = new double[seq_len];
		gamma[i] = new double[seq_len];
	}
	for(i = 0; i < hmm->x_dim; i++)
		for(j = 0; j < hmm->x_dim; j++)
			psi[i][j] = new double[seq_len];


	if(hmm->y_type == CONTINUOUS)
	{
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
			{
				gammagamma[i][j] = new double[seq_len];	
				y_mew_square_exp_vecs[i][j] = new double[seq_len];
			}
	}


	





	double best_model_score = -99999999999;  // the best score to keep ! 
	hmm_model best_hmm; 	


	for(i = 0; i < num_starts; i++)  // try different starting points 
	{
		InitilizeModel(&best_hmm, x_dim, y_dim, y_type, place_flag); // Currently we use random parameters start !!!
	

		// Here maybe use a more clever starting point ... 

		long iter = 0;

		double cur_score = -99999999999; 
		double prev_score = -99999999999.9 - 2*tolerance;

		while(  (iter < max_iters) && (cur_score - prev_score) > tolerance  )
		{
		
			prev_score = cur_score;  // update the previous
				
			// E-step : estimate probabilities	
#ifdef CHECK_TIMING
			clock_t EM_start_time = clock();
#endif
			Compute_b_cond_tabs(&best_hmm, y_vec, y_vec_int, seq_len, b_cond_tabs, y_mew_square_exp_vecs);
#ifdef CHECK_TIMING
			clock_t EM_end_time = clock();
			printf("Compute_b_cond_tabs (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			forward(&best_hmm, seq_len, b_cond_tabs, alpha, scale, scale_cum, scale_exp, &cur_score);
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Forward (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			backward(&best_hmm, seq_len, scale_exp, b_cond_tabs, beta);		
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Backward (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			ComputeMarginalGamma(&best_hmm, y_vec, seq_len, y_mew_square_exp_vecs, alpha, beta, scale_exp, b_cond_tabs, 
				gamma, gammagamma);
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Marginal (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			ComputePsi(&best_hmm, seq_len, alpha, beta, scale, b_cond_tabs, psi);
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Psi (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			// M-step : update model parameters
			UpdateModelParams(&best_hmm, y_vec, y_vec_int, seq_len, y_mew_square_exp_vecs, gamma, gammagamma, psi);	
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Update (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif

			if(iter%10 == 9)
				printf("FINISHED ITER %ld SCORE : %lf  ----\n", iter, cur_score); // PrintModel(&best_hmm); printf("\n---\n Now Realy Finished Iter \n----\n");
			iter++;

		}

		printf("Finished Start %ld Iters %ld cur_score %lf prev_score %lf\n", i, iter, cur_score, prev_score);

		
		if(cur_score > best_model_score)  // here we have improved !!
		{
			best_model_score = cur_score;
			CopyModel(&best_hmm, hmm);
		}

		
		/** Print some ...***
		for(j = 0;j < x_dim; j++)
			for(t=0; t<10; t++) 
				printf("Gam[%ld][%ld] %lf alph= %lf bet= %lf bcon %lf scale %lf\n", 
					j, t, gamma[j][t], alpha[j][t], beta[j][t], b_cond_tabs[j][t], scale[t]);
		**/



	}   // end loop on starting points 

//	printf("Finished  traininingandcopying\n");
	// return all the cums !!! 
	hmm->cum_flag = 1; PermuteModelIncreasingMean(hmm); // compute_model_auxillary_parameters(hmm); 
//	compute_model_auxillary_parameters(hmm); 

	

	// free all allocated buffers
	delete scale;
	delete scale_cum;
	delete scale_exp;
	delete y_vec_int;
	for(i = 0; i < hmm->x_dim; i++)
	{
		delete b_cond_tabs[i];
		delete alpha[i];
		delete beta[i];
		delete gamma[i];
		for(j = 0; j < hmm->x_dim; j++)
			delete psi[i][j];
	}

	if(hmm->y_type == CONTINUOUS)
	{
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
			{
				delete gammagamma[i][j];
				delete y_mew_square_exp_vecs[i][j];
			}
	}


	printf("Finished EM. Best Score is %lf\n", best_model_score);

	return 0; 

}


// Copy all models parameters 
long CopyModel(hmm_model *hmm1, hmm_model *hmm2)
{
	long i, j;

	// copy types 
	hmm2->x_dim = hmm1->x_dim;
	hmm2->y_dim = hmm1->y_dim;
	hmm2->y_type = hmm1->y_type;


	
	// copy parameters
	for(i = 0; i < hmm2->x_dim; i++)
		hmm2->PI[i] = hmm1->PI[i];
	
	for(i = 0; i < hmm2->x_dim; i++)
		for(j = 0; j < hmm2->x_dim; j++)
			hmm2->M[i][j] = hmm1->M[i][j];

	for(i = 0; i < hmm2->x_dim; i++)
		for(j = 0; j < hmm2->y_dim; j++)
			hmm2->N[i][j] = hmm1->N[i][j];

	if(hmm2->y_type == CONTINUOUS)
	{
		for(i = 0; i < hmm2->x_dim; i++)
			for(j = 0; j < hmm2->y_dim; j++)
			{
				hmm2->MEW[i][j] = hmm1->MEW[i][j];
				hmm2->SIGMA[i][j] = hmm1->SIGMA[i][j];
			}
	}


	// complete everything
	hmm2->cum_flag = 1; hmm2->log_flag = 1;
	ComputeModelAuxillaryParameters(hmm2);

	return 0;

}

// Print everything about the model
long PrintModel(hmm_model *hmm)
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
			PrintDoubleVec(hmm->MEW[i], hmm->y_dim);
		printf("Y Std. Matrix :\n");
		for(i = 0; i < hmm->x_dim; i++)
			PrintDoubleVec(hmm->SIGMA[i], hmm->y_dim);
	}



	return 0;

}


// permute the indexes of the model, thus enabling a more fair 
// comparison of the trained and the original model !!!
long PermuteModelIncreasingMean(hmm_model *hmm)
{
	long i, j; 

	double mean_vec[MAX_X_VALS];	
	long indexes[MAX_X_VALS];
	
	double helperx[MAX_X_VALS][MAX_X_VALS];
	double helpery[MAX_X_VALS][MAX_Y_VALS];



	if(hmm->y_type == DISCRETE)  // here permute according to conditional means, assuming the Y's are 0,1,..
	{
		for(i = 0; i < hmm->x_dim; i++)
		{
			mean_vec[i] = 0; 
			for(j = 0; j < hmm->y_dim; j++)
				mean_vec[i] += hmm->N[i][j] * j; 

		}

	}
	else						 // here permute according to the conditional means
	{
		for(i = 0; i < hmm->x_dim; i++)
		{
			mean_vec[i] = 0; 
			for(j = 0; j < hmm->y_dim; j++)
				mean_vec[i] += hmm->N[i][j] * hmm->MEW[i][j];

		}
	}

	
/*	
	printf("start quicksort\n mean : ");
	print_double_vec(mean_vec, hmm->x_dim);
*/
  
	// Now sort according to mean !! 	
	DoQuicksort(mean_vec, hmm->x_dim, indexes);

/*

	printf("ned quicksort\n mean : ");
	print_double_vec(mean_vec, hmm->x_dim);

	printf("\nIndexes :\n");
	print_vec(indexes, hmm->x_dim);
*/

	// change the X's M's variables 
	for(i = 0; i < hmm->x_dim; i++)
		for(j = 0; j < hmm->x_dim; j++)
			helperx[i][j] = hmm->M[indexes[i]][indexes[j]];
	for(i = 0; i < hmm->x_dim; i++)
		for(j = 0; j < hmm->x_dim; j++)
			hmm->M[i][j] = helperx[i][j];

	// Change the Y's N's variables 
	for(i = 0; i < hmm->x_dim; i++)
		for(j = 0; j < hmm->y_dim; j++)
			helpery[i][j] = hmm->N[indexes[i]][j];
	for(i = 0; i < hmm->x_dim; i++)
		for(j = 0; j < hmm->y_dim; j++)
			hmm->N[i][j] = helpery[i][j];

	if(hmm->y_type == CONTINUOUS)  // Change MEW and SIGMA
	{
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				helpery[i][j] = hmm->MEW[indexes[i]][j];
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				hmm->MEW[i][j] = helpery[i][j];

		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				helpery[i][j] = hmm->SIGMA[indexes[i]][j];
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				hmm->SIGMA[i][j] = helpery[i][j];

	}

	// Change everything else 
	ComputeModelAuxillaryParameters(hmm);

	return 0;

}



// Initilize model with given sizes. Unless stated otherwise, parameters are randomized
long InitilizeModel(hmm_model *hmm, long x_dim, long y_dim, long y_type, long place_flag)
{
	long i, j;
	double sumsum;


	// Copy parameters
	hmm->x_dim = x_dim;
	hmm->y_dim = y_dim; 
	hmm->y_type = y_type; 
	hmm->place_flag = place_flag;


//	printf("Start rand ..\n");

	// Start with random parameters unless stated otherwise ..

	
	// Initial distribution on X's
	sumsum = 0; 
	for(i = 0; i < x_dim; i++)
	{
		hmm->PI[i] = randr();
		sumsum += hmm->PI[i];
	}
	for(i = 0; i < x_dim; i++)
		hmm->PI[i] /= sumsum;


	
	// Transition matrices 
	for(i = 0; i < x_dim; i++)
	{
		sumsum = 0; 
		for(j = 0; j < x_dim; j++)
		{
			hmm->M[i][j] = randr();
			sumsum += hmm->M[i][j];
		}
		for(j = 0; j < x_dim; j++)
			hmm->M[i][j] /= sumsum;  // Normalize to make sum equal one !
 
	}

//		printf("Start emission ..\n");

	// Emission matrices
	for(i = 0; i < x_dim; i++)
	{
		sumsum = 0; 
		for(j = 0; j < y_dim; j++)
		{
			hmm->N[i][j] = randr();
			sumsum += hmm->N[i][j];
		}
		for(j = 0; j < y_dim; j++)
			hmm->N[i][j] /= sumsum;  // Normalize to make sum equal one !
 	}
	
	
	// MEW and SIGMA. Currently  MEW ~ U[-1,1] and SIGMA ~ U[0,1]
	if(hmm->y_type == CONTINUOUS)
		for(i = 0; i < x_dim; i++)
			for(j = 0; j < y_dim; j++)
			{
				hmm->MEW[i][j] = 2*randr()-1; 
				hmm->SIGMA[i][j] = randr(); 
			}



	ComputeModelAuxillaryParameters(hmm);


	return 0;
}


// here we only use the model parameters and compute from them the helpfull cumsum and log tables
long ComputeModelAuxillaryParameters(hmm_model *hmm)
{

	long i, j;


	if(hmm->cum_flag)
	{
		// Compute comulative sum for efficiency :
		hmm->PI_cum[0] = hmm->PI[0];
		for(i = 1; i < hmm->x_dim; i++)
			hmm->PI_cum[i] = hmm->PI_cum[i-1] + hmm->PI[i];	

		for(i = 0; i < hmm->x_dim; i++)
		{	

			hmm->M_cum[i][0] = hmm->M[i][0];
			for(j = 1; j < hmm->x_dim; j++)
				hmm->M_cum[i][j] = hmm->M_cum[i][j-1] + hmm->M[i][j];

			hmm->N_cum[i][0] = hmm->N[i][0];
			for(j = 1; j < hmm->y_dim; j++)
				hmm->N_cum[i][j] = hmm->N_cum[i][j-1] + hmm->N[i][j];
		}
	}
	if(hmm->log_flag)
	{

		// Compute logs for efficiency
		for(i = 0; i < hmm->x_dim; i++)
				hmm->PI_log[i] = log(hmm->PI[i]);

		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				hmm->M_log[i][j] = log(hmm->M[i][j]);

		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				hmm->N_log[i][j] = log(hmm->N[i][j]);

		if(hmm->y_type == CONTINUOUS)
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->SIGMA_log[i][j] = log(hmm->SIGMA[i][j]);
	}

	
	for(i = 0; i < hmm->x_dim; i++)
		for(j = 0; j < hmm->y_dim; j++)
			hmm->SIGMA_inv[i][j] = 1.0 / hmm->SIGMA[i][j];

	return 0; 
}



// Give various measures of (dis)similarity between two HMMs
long CompareTwoModels(hmm_model *hmm1, hmm_model *hmm2)
{
	long i, j;


	double M_ssq = 0; double N_ssq = 0; 
	double MEW_ssq = 0, SIGMA_ssq = 0;   // only for continuous
	// First check if two models can be compared at all !!! 
	if(hmm1->y_type  != hmm2->y_type)
	{
		printf("Error !! One HMM is Discrete and one is Continuous. Models cannot be compared !!\n");
		return (-1);
	}

	if(hmm1->y_dim  != hmm2->y_dim)
	{
		printf("Error !! Dimensions of two HMMs don't match !!\n");
		return (-1);
	}


	if(hmm1->x_dim  != hmm2->x_dim)
	{
		printf("Warning !! Dimensions of the X's (hidden) of two HMMs don't match. Avoid some comparisons.\n");
	}
	else
	{
		// Compute sum of square for transition
	
		for(i = 0; i < hmm1->x_dim; i++)
			for(j = 0; j < hmm1->x_dim; j++)
				M_ssq += (hmm1->M[i][j] - hmm2->M[i][j])*(hmm1->M[i][j] - hmm2->M[i][j]);

		
		if(hmm1->y_type  == DISCRETE)
		{
		// Compute sum of square for emission. Currently assume they're discrete 
		
		for(i = 0; i < hmm1->x_dim; i++)
			for(j = 0; j < hmm1->y_dim; j++)
				N_ssq += (hmm1->N[i][j] - hmm2->N[i][j])*(hmm1->N[i][j] - hmm2->N[i][j]);
		}
		else 
		{
			
			for(i = 0; i < hmm1->x_dim; i++)
			{
				for(j = 0; j < hmm1->y_dim; j++)
				{
					MEW_ssq += (hmm1->MEW[i][j] - hmm2->MEW[i][j])*(hmm1->MEW[i][j] - hmm2->MEW[i][j]);   // The means 
					SIGMA_ssq += (hmm1->SIGMA[i][j] - hmm2->SIGMA[i][j])*(hmm1->SIGMA[i][j] - hmm2->SIGMA[i][j]);   // The variances
					N_ssq += (hmm1->N[i][j] - hmm2->N[i][j])*(hmm1->N[i][j] - hmm2->N[i][j]);   // The mixture coefficients
				}
			}
		}
	}


	// Summury of comparison

	printf("Two Models Comparison : \nM_SSQ %lf   N_SSQ %lf\n", M_ssq, N_ssq);
	if(hmm1->y_type  == CONTINUOUS)
		printf("MEW_SSQ %lf   SIGMA_SSQ %lf\n", MEW_ssq, SIGMA_ssq);


	return 0; // everything is OK

}



// Simulate a sequence from a given model. We assume that the vector is already allocated .
// We output also the x vector if someone wants it ...
long SimulateSequenceFromModel(hmm_model *hmm, long seq_len, long *x_vec, double *y_vec)
{

	long j, t; // t denotes time .....

	
	long mix_choose;

//	printf("START SIM\n");

	// first generate the x'th ...
	// Choose the 1st state according to the initial distribution 
	double randy = randr();	
	for(j = 0; j < hmm->x_dim; j++)
		if(randy < hmm->PI_cum[j])
		{
			x_vec[0] = j;
			break;
		}


	// Choose all other states
	for(t = 1; t < seq_len; t++)
	{


		randy = randr();
		for(j = 0; j < hmm->x_dim; j++)
		if(randy < hmm->M_cum[x_vec[t-1]][j])
		{
			x_vec[t] = j;
			break;
		}
	}

//	printf("FINISHED X's\n");


	// Now simulate the y's
	if(hmm->y_type == DISCRETE)	
		for(t = 0; t < seq_len; t++)
		{
			randy = randr();
			for(j = 0; j < hmm->y_dim; j++)
				if(randy < hmm->N_cum[x_vec[t]][j])
				{
					y_vec[t] = j;
					break;
				}
		}
	else  // Here Mixture of Gaussians
		for(t = 0; t < seq_len; t++)
		{
			// First randomize which of the mixtures to take
			randy = randr();
			for(j = 0; j < hmm->y_dim; j++)
				if(randy < hmm->N_cum[x_vec[t]][j])
				{
					mix_choose = j;
					break;
				}

			// Now randomize according to a gaussian distribution 
			y_vec[t] = Gaussian(hmm->MEW[x_vec[t]][mix_choose], hmm->SIGMA[x_vec[t]][mix_choose]);

		}


//	printf("FINISHED Y's\n");

	return 0;

}





#endif