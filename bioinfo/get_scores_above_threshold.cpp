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

  long is_double; 
	
  double **scores;
  float **scores_single;
  long num_seqs; 
  long seqs_len; 
  double threshold;
  float threshold_single;
  long *out_x, *out_y;
  double *out_scores;
  float *out_scores_single;

  double *in, *inn, *out;
  float *in_single, *out_single; 
  long *seqs_lens; // If we have many sequences with different lengthes   


 // printf("Start get_scores_above_threshold\n");


  /* Check for proper number of arguments. Should be Five 5 (last one is double flag)*/
  if(nrhs != 5) {
    mexErrMsgTxt("Usage: [dat] [K]");
  }  

  in = mxGetPr(prhs[4]);   // Get the is_double flag
  is_double	= long(in[0]);	  


  if(is_double)
  {
     in = mxGetPr(prhs[0]);   // Get the scores matrix
     num_seqs = mxGetM(prhs[0]);    // Get the number of sequences (genes)
     seqs_len = mxGetN(prhs[0]);    // Get the length of the sequences (assumed to be the same for all)	  

     scores = new double*[num_seqs];
     for(i = 0; i < num_seqs; i++)
	    scores[i] = new double[seqs_len];  // allocate memory for the scores
     for(i = 0; i < num_seqs; i++)
	    for(j = 0; j < seqs_len; j++)
		   scores[i][j] = in[j *  num_seqs + i  /*i*num_seqs+j*/];  // Copy the sequences matrix (Note that they should be packed for now !!!)
  }
  else
  {
	 in_single = (float *)(mxGetData(prhs[0]));		// Get the scores matrix (singles)
     num_seqs = mxGetM(prhs[0]);    // Get the number of sequences (genes)
     seqs_len = mxGetN(prhs[0]);    // Get the length of the sequences (assumed to be the same for all)	  

     scores_single = new float*[num_seqs];
     for(i = 0; i < num_seqs; i++)
	    scores_single[i] = new float[seqs_len];  // allocate memory for the scores
     for(i = 0; i < num_seqs; i++)
	    for(j = 0; j < seqs_len; j++)
		   scores_single[i][j] = in_single[j *  num_seqs + i  /*i*num_seqs+j*/];  // Copy the sequences matrix (Note that they should be packed for now !!!)
  }

  seqs_lens = new long[num_seqs]; // allocate the lengths array
  for(i = 0; i < num_seqs; i++)
	  seqs_lens[i] = seqs_len;    // All genes are assumed to have the same length

  in = mxGetPr(prhs[1]);   // Get the threshold   	  
  threshold	= in[0];	   // to be used 
  
  in = mxGetPr(prhs[2]);   // Get the max overlap distance   	  
  long remove_overlap_distance = long(in[0]);	   // to be used 
 
  in = mxGetPr(prhs[3]);   // Get the smart_thresh flag   	  
  long smart_thresh = long(in[0]);	   // to be used 

  // check what is good ind	
  long good_ind = 0;


	
  // First go over all scores to see how many are above threshold
  for(g = 0; g < num_seqs; g++)
	 for(j = 0; j < seqs_lens[g]; j++)
		if(scores[g][j] >= threshold)
		   good_ind++;

  // allocate enough memory for outputs of the function 
  out_x = new long[good_ind];
  out_y = new long[good_ind];
  if(is_double)
     out_scores = new double[good_ind];
  else
	  out_scores_single = new float[good_ind];

  good_ind = get_scores_above_threshold(is_double, scores, scores_single, num_seqs, seqs_lens, threshold, threshold_single,
										remove_overlap_distance, smart_thresh, 
										out_x, out_y, out_scores, out_scores_single); // Those are allocated OUTSIDE the function

  /**
  long get_scores_above_threshold(long is_double, double *scores[MAX_NUM_SEQS], float *scores_single[MAX_NUM_SEQS], 
								long num_seqs, long *seqs_lens, double threshold, float threshold_single,
								long remove_overlap_distance, long smart_thresh, 
								long *out_x, long *out_y, 
								double *out_scores, float *out_scores_single); // Those are allocated OUTSIDE the function
**/


  // printf("finished call get_scores_above_threshold555\n");

	// The three outputs of the function 
	plhs[0] = mxCreateDoubleMatrix(good_ind, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(good_ind, 1, mxREAL);
	


  out = mxGetPr(plhs[0]); // output the x's array 
  for(i=0; i<good_ind; i++)
	  out[i] = out_x[i];
  out = mxGetPr(plhs[1]); // output the y's array 
  for(i=0; i<good_ind; i++)
	  out[i] = out_y[i];
  if(is_double)
  {
     plhs[2] = mxCreateDoubleMatrix(good_ind, 1, mxREAL);
     out = mxGetPr(plhs[2]); // output the scores's array 
     for(i=0; i<good_ind; i++)
	    out[i] = out_scores[i];
  }
  else
  {
	 plhs[2] = mxCreateNumericMatrix(good_ind, 1, mxSINGLE_CLASS, mxREAL);
	  out_single = (float *)(mxGetData(plhs[2]));  // output the scores's array  (singles)
	  for(i=0; i<good_ind; i++)
	     out_single[i] = out_scores_single[i];
  }

  // free memory
  if(is_double)
  {
     for(i=0; i<num_seqs; i++) 
	    delete scores[i];
     delete scores;
  }
  else
  {
     for(i=0; i<num_seqs; i++) 
	    delete scores_single[i];
     delete scores_single;
  }
  delete seqs_lens;
  delete out_x;
  delete out_y;
  delete out_scores;


}
//#endif // MEX_TO_MATLAB



/////////////////
// A list of functions dealing with DNA
// Note : Here we deal with DNA in COMPRESSED form. 
// This will be the standard from now on !!!!!!!
// We will use chars of 8 bits. The 2 lsb will only be considered
// Each 2 bits will be one DNA base. 
// This will be : (from Now on !!!!)
// A - 00
// T - 01
// C - 10
// G - 11
// a DNA sequence will always be an array (of chars ???), with the first char saying the length !!!
// There should always be enough place allocated
// This causes a 4-times extra memory, but what the f***
/////////////////





/**************************************************************************/

/******************************
 * qsort.c                    *
 * quicksort algorithm in C   *
 * Daniel J. Schultz          *
 ******************************/

// inner split routine
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










// Get an array of scores, and output the locations & scores of the ones which have past the threshold.
// In the future it is possible to let the routine select it's own threshold in some smart way.
long get_scores_above_threshold(double *scores[MAX_NUM_SEQS], long num_seqs, long *seqs_lens, double threshold, 
								long remove_overlap_distance, long smart_thresh, 
								long *out_x, long *out_y, double *out_scores) // Those are allocated OUTSIDE the function
{
	
	
	long i, j, k, g;
	
	
	long dist = remove_overlap_distance; 

	long count_above = 0;
	long counter = 0;
	double local_threshold;

	local_threshold = threshold;
//	long *indexes; 

//	printf("STARTING get_scores_above_threshold\n");


	// Now determine the threshold if needed. Here at the beginning threshold is the FRACTION !!! 
//		if(smart_threshold == 1)
//		for(g = 0; g < num_seqs; g++)
//		{
//			DoQuicksort(scores[g], seqs_lens[g], indexes)
//		}

			
		// First go over all scores to see how many are above threshold
		for(g = 0; g < num_seqs; g++)
			for(j = 0; j < seqs_lens[g]; j++)
				if(scores[g][j] >= local_threshold)
					count_above++;



//	printf("DONE counting\n");


			
/*** We assume that those arrays are already allocated from outside ..
	out_x = new long[count_above];
	out_y = new long[count_above];
	out_scores = new double[count_above];
***/

			///printf("Number of genes here : %ld , Gene's Length here : %ld\n", num_seqs, seqs_lens[0]);

		
	// Now loop again over all scores to fill the new arrays
	for(g = 0; g < num_seqs; g++)
		for(j = 0; j < seqs_lens[g]; j++)
			if(scores[g][j] >= local_threshold)
			{
//				printf("g+1   %ld    j+1 %ld  --- ", g+1, j+1);
				out_x[counter] = g+1; 
				out_y[counter] = j+1;   // we do plus 1 because we go to matlab ...
//				printf("x   %ld    y %ld   NOW \n", out_x[counter], out_y[counter]);
				out_scores[counter++] = scores[g][j];
//				printf("x   %ld    y %ld   score %lf\n", out_x[counter-1], out_y[counter-1], out_scores[counter-1]);
			}
//	printf("DONE counting againnnn\n");

			
	long count_good = count_above; // The number of good scores left

///	dist = -1; // try not to remove overlaps at first ...
	if(dist != (-1))
	{
		long *good_ind = new long[count_above];
		for(i = 0; i < count_above; i++)
			good_ind[i] = 1;

		for(i = 0; i < count_above; i++)
		{
			for(k = i+1; k < count_above; k++)  // go forward
			{
				if(out_x[k] != out_x[i])  // we got to some other gene
					break; 
				if(out_y[k] > out_y[i] + dist)  // here we are more than dist away
					break;
				if(out_scores[k] > out_scores[i])  // here we got a better scores !!!
				{
					good_ind[i] = 0; count_good--;
					break;
				}
			}


			if(good_ind[i])
				for(k = i-1; k >=0; k--)  // go also backwards if needed
				{
					if(out_x[k] != out_x[i])  // we got to some other gene
						break; 
					if(out_y[k] < out_y[i] - dist)  // here we are more than dist away
						break;
					if(out_scores[k] > out_scores[i])  // here we got a better scores !!!
					{
						good_ind[i] = 0; count_good--;
						break;
					}
				}

		}


		// Now that we know the good indices, we can set the array .
		long *new_x = new long[count_good];
		long *new_y = new long[count_good];
		double *new_scores = new double[count_good];

		counter = 0;
		for(i = 0; i < count_above; i++)
			if(good_ind[i])
			{
				new_x[counter] = out_x[i];
				new_y[counter] = out_y[i];
				new_scores[counter++] = out_scores[i];
			}


		// Set the pointers back  - We cannot set pointers since we will 'lose' them, so we copy whole the values
/***		out_x = new_x;
		out_y = new_y;
		out_scores = new_scores;
		***/
		for(i = 0; i < count_good; i++)
		{
			out_x[i] = new_x[i];
			out_y[i] = new_y[i];
			out_scores[i] = new_scores[i];
		}

		delete good_ind;
		delete new_x;
		delete new_y;
		delete new_scores;
	}





//	printf("Returned %ld GOOD SITES \n", count_good);
	return count_good;


}



