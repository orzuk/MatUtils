#include <stdlib.h>
#include <stdio.h>
//#include <iostream.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "general.h"
#include "dna_utils_new.h"
#include "markov.h"


///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For Yuval's Project - 11.11.2003
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int nRows, nCols, i, j, K, g;
	
  double **scores; 
  float **scores_single; // we choose later if to use singles or doubles 
  long num_seqs; 
  long seqs_len; 
  double threshold; float threshold_single;
  long *out_x, *out_y;
  double *out_scores;  float *out_scores_single; 

  double *in, *inn, *out;
  float *in_single, *out_single;
  long *seqs_lens; // If we have many sequences with different lengthes   

  /* Check for proper number of arguments. Should be Six 6 */
  if(nrhs != 6) {
    mexErrMsgTxt("Usage: [dat] [K]");
  }
  
  in = mxGetPr(prhs[5]);
  long is_double = long(in[0]); // New! flag saying if we read a single or double array 


//  seqs_len = mxGetN(prhs[0]);    // Get the length of the sequences (assumed to be the same for all)	  

  in = mxGetPr(prhs[1]);   // Get the scores lengths
  num_seqs = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));    // Get the number of sequences (genes)
  seqs_lens = new long[num_seqs]; // allocate the lengths array
  for(i = 0; i < num_seqs; i++)
	  seqs_lens[i] = long(in[i]);    // Read genes lengths from input
//   printf("num seqs is %ld IS DOUBLE START IS: %ld\n", num_seqs, is_double);

  // Now, that we know the lengths, we can copy the scores 
  long ind=0;
  if(is_double)
  {  
	  in = mxGetPr(prhs[0]);   // Get the scores matrix (doubles)
	  scores = new double*[num_seqs];
	  for(i = 0; i < num_seqs; i++)
		scores[i] = new double[seqs_lens[i]];  // allocate memory for the scores doubles
	  for(i = 0; i < num_seqs; i++)
	  {
		for(j = 0; j < seqs_lens[i]; j++)
			scores[i][j] = in[ind+j];  // Copy the sequences matrix (Note that they should be packed for now !!!)
		ind += seqs_lens[i];
	  }
	  in = mxGetPr(prhs[2]);   // Get the threshold to be used  	  
	  threshold	= in[0];	   
  }
  else
  {
	  in_single = (float *)(mxGetData(prhs[0]));		// Get the scores matrix (singles)
	  scores_single = new float*[num_seqs];
	  for(i = 0; i < num_seqs; i++)
	     scores_single[i] = new float[seqs_lens[i]];  // allocate memory for the scores singles
	  for(i = 0; i < num_seqs; i++)
	  {
		for(j = 0; j < seqs_lens[i]; j++)
			scores_single[i][j] = in_single[ind+j];  // Copy the sequences matrix (Note that they should be packed for now !!!)
		ind += seqs_lens[i];
	  }
	  in_single = (float *)(mxGetData(prhs[2]));	 // Get the threshold to be used  	  
	  threshold_single = in_single[0];	   
  }	  

  in = mxGetPr(prhs[3]);   // Get the max overlap distance   	  
  long remove_overlap_distance = long(in[0]);	   // to be used 
  in = mxGetPr(prhs[4]);   // Get the smart_thresh flag   	  
  long smart_thresh = long(in[0]);	   // to be used 

//  printf("Remove Overlap Distance Read: %ld\n", remove_overlap_distance);

  // check what is good ind	
  long good_ind = 0;
	
  // First go over all scores to see how many are above threshold
/////  if(smart_thresh == 0) 
  long total_len = 0;

  if(is_double)
  {
	for(g = 0; g < num_seqs; g++)
		for(j = 0; j < seqs_lens[g]; j++)
		{
			if(scores[g][j] >= threshold)
				good_ind++;
			total_len++;
		}
  }
  else
  {
	for(g = 0; g < num_seqs; g++)
		for(j = 0; j < seqs_lens[g]; j++)
		{
			if(scores_single[g][j] >= threshold_single)
				good_ind++;
			total_len++;
		}
  }

//  printf("good ind is %ld threshold_single is %lf\n", 	good_ind, threshold_single); 

//	printf("NUMBER OF SEQS %ld TOTAL NUMBER OF CANDIDATE SCORES : %ld  C Threshold %lf Init Good Ind %ld\n",
//		num_seqs, total_len, threshold, good_ind);
//	fflush(stdout);

////  else  // here the threshold is given as a fraction !!! 
////	  for(g = 0; g < num_seqs; g++)
////		  good_ind = good_index + ceil(seqs_lens[g]*threshold); 

  // allocate enough memory for outputs of the function 
  out_x = new long[good_ind];
  out_y = new long[good_ind];
  if(is_double)
	out_scores = new double[good_ind];
  else
	out_scores_single = new float[good_ind];


/*****
	long get_scores_above_threshold(double *scores[MAX_NUM_SEQS], long num_seqs, long *seqs_lens, double threshold, 
								long remove_overlap_distance, long smart_thresh, 
								long *out_x, long *out_y, double *out_scores) // Those are allocated OUTSIDE the function
*****/
///	return; // EXITING 

 //   printf("call get_scores_above_threshold\n");
//    for(i = 0; i < num_seqs; i++)
//	  for(j = 0; j < seqs_lens[i]; j++)
//		  printf("Before Func: Read element %ld and its %lf \n", ind+j, scores[i][j]);
//    printf("Threshold is: %lf   Overlap dist is: %ld \n", threshold, remove_overlap_distance);	
/***/
	good_ind = get_scores_above_threshold(is_double, scores, scores_single, num_seqs, seqs_lens, threshold, threshold_single,
									remove_overlap_distance, smart_thresh, 
									out_x, out_y, out_scores, out_scores_single); // Those are allocated OUTSIDE the function
/***/

   if(is_double) // all outputs are double 
   {
	  plhs[0] = mxCreateDoubleMatrix(good_ind, 1, mxREAL); 	// The three outputs of the function 
	  plhs[1] = mxCreateDoubleMatrix(good_ind, 1, mxREAL);
	  plhs[2] = mxCreateDoubleMatrix(good_ind, 1, mxREAL);

	  out = mxGetPr(plhs[0]); // output the x's array 
	  for(i=0; i<good_ind; i++)
	     out[i] = out_x[i];
	  out = mxGetPr(plhs[1]); // output the y's array 
	  for(i=0; i<good_ind; i++)
	     out[i] = out_y[i];
	  out = mxGetPr(plhs[2]); // output the scores's array 
	  for(i=0; i<good_ind; i++)
	     out[i] = out_scores[i];
   }
   else // here we work with singles 
   {
      plhs[0] = mxCreateDoubleMatrix(good_ind, 1, mxREAL); 	// The three outputs of the function 
	  plhs[1] = mxCreateDoubleMatrix(good_ind, 1, mxREAL);
	  plhs[2] = mxCreateNumericMatrix(good_ind, 1, mxSINGLE_CLASS, mxREAL);

	  out = mxGetPr(plhs[0]); // output the x's array 
	  for(i=0; i<good_ind; i++)
	    out[i] = out_x[i];
	  out = mxGetPr(plhs[1]); // output the y's array 
	  for(i=0; i<good_ind; i++)
	    out[i] = out_y[i];
	  out_single = (float *)(mxGetData(plhs[2]));  // output the scores's array  (singles)
	  for(i=0; i<good_ind; i++)
	     out_single[i] = out_scores_single[i];
   }


  // free memory
  delete seqs_lens;
  delete out_x;
  delete out_y;
  if(is_double)
  {
	  for(i=0; i<num_seqs; i++) 
		delete scores[i];
	  delete scores;
	  delete out_scores;
  }
  else
  {
	  for(i=0; i<num_seqs; i++) 
		delete scores_single[i];
	  delete scores_single;
	  delete out_scores_single;
  }
}
//#endif // MEX_TO_MATLAB




/**************************************************************************/

/******************************
 * qsort.c                    *
 * quicksort algorithm in C   *
 * Daniel J. Schultz          *
 ******************************/

// inner split routine
/************
long split(double *vals, long first, long last, long *indexes);
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
long quicksort(double *vals, long first, long last, long *indexes);
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
long DoQuicksort(double *vals, long len, long *indexes);
long DoQuicksort(double *vals, long len, long *indexes)
{
	long i;

	for(i = 0; i < len; i++)
		indexes[i] = i;  // init the indexes 
		
	quicksort(vals, 0, len-1, indexes);

	return 0; 
}
**************/










