#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "../hmp/general.h"
#include "../hmp/hmm_chrom_funcs.h"


long PrintModel2(hmm_model *hmm);
long PrintDoubleVec2(double *vec, long len);

///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For HMM Chromosome project - 01.08.2004
/// Fit a MoG model to data
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
	int nRows, nCols, i, j, K, t, max_iters, num_starts;

	long num_points; 
	long x_dim, y_dim, place_flag, miss_data, num_iters;
	double tolerance;
  
	double *in, *out, *data;


	double p[MAX_NUM_GAUSSIANS], m[MAX_NUM_GAUSSIANS], s[MAX_NUM_GAUSSIANS];
	double init_p[MAX_NUM_GAUSSIANS], init_m[MAX_NUM_GAUSSIANS], init_s[MAX_NUM_GAUSSIANS];
	double log_like;

	 /* Check for proper number of arguments. Should be six 6 */
	if(nrhs != 6) 
		mexErrMsgTxt("Usage: [dat] [K]");


  /* Parameters should be : 
	1. Data vector
	2. Number of gaussians 
	3. Number of iterations to perform
	4. Initial P initial probabilities vector
	5. Initial M mean vector
	6. Initial S standard error matrix

 ---   Note : Dimensions are already IN the model !!! ---
	*/

///////////////////////////////////////////////////////////////
//// Read the data ...
///////////////////////////////////////////////////////////////
	in = mxGetPr(prhs[0]);		// Get the data vector
	num_points = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the data length 
	data = new double[num_points];
	for(i = 0; i < num_points; i++)
		data[i] = in[i];
	in = mxGetPr(prhs[1]);		// Get the number of gaussians
	long num_gaussians = long(in[0]);	// Get the number of gaussians

	in = mxGetPr(prhs[2]);		// Get the number of iterations
	num_iters = long(in[0]);	// Get the number of iterations
	in = mxGetPr(prhs[3]);		// Get the init probs vector
	for(i = 0; i < num_gaussians; i++)
		init_p[i] = in[i];      // Get the init probs vector
	in = mxGetPr(prhs[4]);		// Get the init mu vector
	for(i = 0; i < num_gaussians; i++)
		init_m[i] = in[i];      // Get the init mu vector
	in = mxGetPr(prhs[5]);		// Get the init sigma vector
	for(i = 0; i < num_gaussians; i++)
		init_s[i] = in[i];      // Get the init sigma vector


//	printf("num_points %ld num_gaussians %ld num_iters %ld\n", 
//		long(num_points), long(num_gaussians), num_iters);
	// Run the MoG function
	double mog_ret = MixtureOfGaussiansGivenInit(data, long(num_points), long(num_gaussians), num_iters, 
		 init_p, init_m, init_s, 
		 p, m, s, &log_like);
//	printf("Finished MOG, LogLike = %lf\n", log_like);

	// return the output values
	plhs[0] = mxCreateDoubleMatrix(1, num_gaussians, mxREAL);  // The probs vector 
	out = mxGetPr(plhs[0]); // output probs vec
	for(i = 0; i < num_gaussians; i++)
		out[i] = p[i]; 
	plhs[1] = mxCreateDoubleMatrix(1, num_gaussians, mxREAL);  // The means vector 
	out = mxGetPr(plhs[1]); // output means vec
	for(i = 0; i < num_gaussians; i++)
		out[i] = m[i]; 
	plhs[2] = mxCreateDoubleMatrix(1, num_gaussians, mxREAL);  // The stds vector 
	out = mxGetPr(plhs[2]); // output stds vec
	for(i = 0; i < num_gaussians; i++)
		out[i] = s[i]; 
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);  // The log_like vector 
	out = mxGetPr(plhs[3]); // output log_like vec
	out[0] = log_like; 



	delete data; // delete memory
}
//#endif // MEX_TO_MATLAB
