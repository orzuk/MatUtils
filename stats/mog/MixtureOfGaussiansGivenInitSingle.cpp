#include <stdlib.h>
#include <stdio.h>
// #include <iostream.h>
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
/// Fit a MoG model to data
/// The same as MixtureOfGaussianGivenInit but for Single data input
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
	long nRows, nCols, i, j, K, t, max_iters, num_starts;

	long num_points; 
	long x_dim, y_dim, place_flag, miss_data, num_iters;
	float tolerance;
  
	float *in, *out, *data;


	float p[MAX_NUM_GAUSSIANS], m[MAX_NUM_GAUSSIANS], s[MAX_NUM_GAUSSIANS];
	float init_p[MAX_NUM_GAUSSIANS], init_m[MAX_NUM_GAUSSIANS], init_s[MAX_NUM_GAUSSIANS];
	float log_like;

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
	in = (float *)(mxGetData(prhs[0]));		// Get the data vector   // mxGetPr
	num_points = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));     // Get the data length 
	data = new float[num_points];
//	double x=0.0; 	printf("sum x is %lf\n", x);
	for(i = 0; i < num_points; i++)
		data[i] = in[i];

	in = (float *)(mxGetData(prhs[1]));		// Get the number of gaussians
	long num_gaussians = long(in[0]);	// Get the number of gaussians

	in = (float *)(mxGetData(prhs[2]));		// Get the number of iterations
	num_iters = long(in[0]);	// Get the number of iterations
	in = (float *)(mxGetData(prhs[3]));		// Get the init probs vector
	for(i = 0; i < num_gaussians; i++)
		init_p[i] = in[i];      // Get the init probs vector
	in = (float *)(mxGetData(prhs[4]));		// Get the init mu vector
	for(i = 0; i < num_gaussians; i++)
		init_m[i] = in[i];      // Get the init mu vector
	in = (float *)(mxGetData(prhs[5]));		// Get the init sigma vector
	for(i = 0; i < num_gaussians; i++)
		init_s[i] = in[i];      // Get the init sigma vector


	// Run the MoG function
	float mog_ret = MixtureOfGaussiansGivenInitSingle(data, long(num_points), long(num_gaussians), num_iters, 
		 init_p, init_m, init_s, 
		 p, m, s, &log_like);  

	printf("Finished MOG, LogLike = %lf\n", log_like);
	for(i = 0; i < num_gaussians; i++)
		printf("Gaussian: %ld  Prob: %lf Mean: %lf  Std: %lf\n", i, p[i], m[i], s[i]); 

	// return the output values
	plhs[0] = mxCreateNumericMatrix(1, num_gaussians, mxSINGLE_CLASS, mxREAL);  // The probs vector 
	out = (float *)(mxGetData(plhs[0])); // output probs vec
	for(i = 0; i < num_gaussians; i++)
		out[i] = p[i]; 
	plhs[1] = mxCreateNumericMatrix(1, num_gaussians, mxSINGLE_CLASS, mxREAL);  // The means vector 
	out = (float *)(mxGetData(plhs[1])); // output means vec
	for(i = 0; i < num_gaussians; i++)
		out[i] = m[i]; 
	plhs[2] = mxCreateNumericMatrix(1, num_gaussians, mxSINGLE_CLASS, mxREAL);  // The stds vector 
	out = (float *)(mxGetData(plhs[2])); // output stds vec
	for(i = 0; i < num_gaussians; i++)
		out[i] = s[i]; 
	plhs[3] = mxCreateNumericMatrix(1, 1, mxSINGLE_CLASS, mxREAL);  // The log_like vector 
	out = (float *)(mxGetData(plhs[3])); // output log_like vec
	out[0] = log_like; 




	delete data; // delete memory
}
//#endif // MEX_TO_MATLAB
