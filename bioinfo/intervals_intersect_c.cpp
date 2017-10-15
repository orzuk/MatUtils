#include <stdlib.h>
#include <stdio.h>
//#include <iostream.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "general.h"
#include "dna_utils_new.h"
#include "hmm_chrom_funcs.h"
#include "markov.h"


///////////////////////////////////////////////////////////////////////////////////////
/// Intersect two sets of intervals - 07.10.2009
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
	long i, j, n1, n2;	
	long inter_ctr, is_sorted;
	double *a, *b;
	double *in, *out;
	double *start_pos1, *end_pos1, *start_pos2, *end_pos2;

	/* Check for proper number of arguments. Should be Five 5 */
	if(nrhs != 5) {
		mexErrMsgTxt("Usage: [dat] [K]");
	}  

	n1 = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));    // Get the number of sequences (genes)
//	printf("Allocating 1st Input: %ld \n", n1); 
	
	start_pos1 = (double *) mxMalloc( n1 * sizeof(double));
	end_pos1 = (double *) mxMalloc( n1 * sizeof(double));
/*	start_pos1 = new double[n1];
	end_pos1 = new double[n1];
	*/
	in = mxGetPr(prhs[0]);   // Get the first values vector matrix
	for(i = 0; i < n1; i++)
		start_pos1[i] = in[i];  // copy first start values
	in = mxGetPr(prhs[1]);   // Get the first values vector matrix
	for(i = 0; i < n1; i++)
		end_pos1[i] = in[i];  // copy first end values

	n2 = MAX(mxGetM(prhs[2]), mxGetN(prhs[2]));    // Get the number of sequences (genes)
//	printf("Allocating 2nd Input: %ld \n", n2); 

	start_pos2 = (double *) mxMalloc( n2 * sizeof(double));
	end_pos2 = (double *) mxMalloc( n2 * sizeof(double));
/*
	start_pos2 = new double[n2];
	end_pos2 = new double[n2];
*/
	in = mxGetPr(prhs[2]);   // Get the first values vector matrix
	for(i = 0; i < n2; i++)
		start_pos2[i] = in[i];  // copy first start values
	in = mxGetPr(prhs[3]);   // Get the first values vector matrix
	for(i = 0; i < n2; i++)
		end_pos2[i] = in[i];  // copy first end values
	in = mxGetPr(prhs[4]);   // Get the sort flag
	is_sorted = in[0]; 

	// print  values
/*	printf("n1 is %ld , n2 is %ld\n", n1, n2);
	for(i=0; i<n1; i++)
		printf("read edges1: (%ld, %ld)\n", start_pos1[i], end_pos1[i]); 
	for(i=0; i<n2; i++)
		printf("read edges2: (%ld, %ld)\n", start_pos2[i], end_pos2[i]); 
*/

//	printf("Allocating Intersect Output: %ld \n", n1+n2); 
	long max_len = MAX(n1,n2); // (double *) mxMalloc( d * sizeof(double));

	double *intersect_start_pos = (double *) mxMalloc( (n1+n2) * sizeof(double));
	double *intersect_end_pos = (double *) mxMalloc( (n1+n2) * sizeof(double));
	double *intersect_inds1 = (double *) mxMalloc( (n1+n2) * sizeof(double));
	double *intersect_inds2 = (double *) mxMalloc( (n1+n2) * sizeof(double));

/*
	double *intersect_start_pos = new double[n1+n2];
	double *intersect_end_pos = new double[n1+n2];
	double *intersect_inds1 = new double[n1+n2];
	double *intersect_inds2 = new double[n1+n2];
*/
	is_sorted = 1; 
	intervals_intersect(intersect_start_pos, intersect_end_pos, intersect_inds1, intersect_inds2, &inter_ctr, 
		start_pos1, end_pos1, start_pos2, end_pos2, n1, n2, is_sorted); // call function

	mxFree(start_pos1); mxFree(start_pos2); mxFree(end_pos1); mxFree(end_pos2);
//	delete start_pos1, start_pos2, end_pos1, end_pos2; // free more memory


	// The two outputs of the function 
	plhs[0] = mxCreateDoubleMatrix(inter_ctr, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(inter_ctr, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(inter_ctr, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(inter_ctr, 1, mxREAL);

//	printf("Output inter_ctr is %ld \n", inter_ctr);
	
	out = mxGetPr(plhs[0]); // output the starts array 
	for(i=0; i<inter_ctr; i++)
		out[i] = intersect_start_pos[i];
	out = mxGetPr(plhs[1]); // output the ends array 
	for(i=0; i<inter_ctr; i++)
		out[i] = intersect_end_pos[i];
	out = mxGetPr(plhs[2]); // output the inds1 array 
	for(i=0; i<inter_ctr; i++)
		out[i] = intersect_inds1[i];
	out = mxGetPr(plhs[3]); // output the inds2 array 
	for(i=0; i<inter_ctr; i++)
		out[i] = intersect_inds2[i];


/*	mxSetPr(plhs[0], intersect_start_pos); 
	mxSetPr(plhs[1], intersect_end_pos); 
	mxSetPr(plhs[2], intersect_inds1); 
	mxSetPr(plhs[3], intersect_inds2); 
*/

	mxFree(intersect_start_pos); mxFree(intersect_end_pos); mxFree(intersect_inds1); mxFree(intersect_inds2); 	// free memory
//	delete intersect_start_pos, intersect_end_pos, intersect_inds1, intersect_inds2; 	// free memory
//	delete out, plhs; 
}
//#endif // MEX_TO_MATLAB


// Find for each element of a the closest element of b (provide speedup for matlab function)
// We assume that they are already sorted and allocated 

/*** Matlab code here: 
n1 = length(start_pos1); n2 = length(start_pos2);
tmp_max = max(start_pos1, end_pos1); tmp_min = min(start_pos1, end_pos1); start_pos1 = tmp_min; end_pos1 = tmp_max;
tmp_max = max(start_pos2, end_pos2); tmp_min = min(start_pos2, end_pos2); start_pos2 = tmp_min; end_pos2 = tmp_max; % make sure end > start
[start_pos1 sort_perm1] = sort(start_pos1); end_pos1 = end_pos1(sort_perm1); % sort the first positions vec by start
[start_pos2 sort_perm2] = sort(start_pos2); end_pos2 = end_pos2(sort_perm2); % sort the second positions vec by start

prev_ctr2 = 1; % counters following both regions

intersect_start_pos = zeros(1, min(n1,n2)); intersect_end_pos = zeros(1, min(n1,n2)); inter_ctr = 1;
intersect_inds1 = zeros(1, min(n1,n2)); intersect_inds2 = zeros(1, min(n1,n2));
for ctr1 = 1:n1
    ctr2 = find(start_pos1(ctr1) <= end_pos2(prev_ctr2:end), 1); % find the first chance for intersection
    if(~isempty(ctr2))
        inter_start = max(start_pos1(ctr1), start_pos2(prev_ctr2+ctr2-1));
        inter_end = min(end_pos1(ctr1), end_pos2(prev_ctr2+ctr2-1));
        while(inter_start <= inter_end)
            intersect_start_pos(inter_ctr) = inter_start;
            intersect_end_pos(inter_ctr) = inter_end;
            intersect_inds1(inter_ctr) = ctr1;   % New: get also the indexes
            intersect_inds2(inter_ctr) = prev_ctr2+ctr2-1;
            inter_ctr = inter_ctr+1;
            if(prev_ctr2+ctr2 <= n2) % make sure we didn't reach the end of second array
                inter_start = max(start_pos1(ctr1), start_pos2(prev_ctr2+ctr2));
                inter_end = min(end_pos1(ctr1), end_pos2(prev_ctr2+ctr2));
            else
                inter_start = 1; inter_end = -1;
            end
            ctr2=ctr2+1;
        end
        prev_ctr2 = max(prev_ctr2, prev_ctr2+ctr2-2); % shift previous counter to the right
    end
end
intersect_start_pos = intersect_start_pos(1:inter_ctr-1); intersect_end_pos = intersect_end_pos(1:inter_ctr-1); % take only relevant inds
intersect_inds1 = intersect_inds1(1:inter_ctr-1); intersect_inds2 = intersect_inds2(1:inter_ctr-1); % take only relevant inds
intersect_inds1 = sort_perm1(intersect_inds1); intersect_inds2 = sort_perm2(intersect_inds2);  % transfer back to original indices
***/
long intervals_intersect(double *intersect_start_pos, double *intersect_end_pos, double *intersect_inds1, double *intersect_inds2, long *inter_ctr,
		double *start_pos1, double *end_pos1, double *start_pos2, double *end_pos2, long n1, long n2, long is_sorted) // call function
{
	long i, j, k, ctr1, ctr2;
	long max_len = MAX(n1,n2);
	double tmp_max, tmp_min; 
	long *sort_perm1;
	long *sort_perm2;
	double *tmp_pos;


	if(!is_sorted) // avoid sort when not needed 
	{
		printf("Allocating Sort Perm: %ld %ld %ld\n", n1,n2,n1+n2); 
		sort_perm1 = (long *) mxMalloc( n1 * sizeof(long));
		sort_perm2 = (long *) mxMalloc( n2 * sizeof(long));
		tmp_pos = (double *) mxMalloc( (n1+n2) * sizeof(double));
		/*
		sort_perm1 = new long[n1];
		sort_perm2 = new long[n2];
		tmp_pos = new double[n1+n2];
		*/
		for(i=0; i<n1; i++)  // make sure end > start
		{
			tmp_max = MAX(start_pos1[i], end_pos1[i]);
			tmp_min = MIN(start_pos1[i], end_pos1[i]);
			start_pos1[i] = tmp_min; end_pos1[i] = tmp_max;
		}	
		for(i=0; i<n2; i++)
		{
			tmp_max = MAX(start_pos2[i], end_pos2[i]);
			tmp_min = MIN(start_pos2[i], end_pos2[i]);
			start_pos2[i] = tmp_min; end_pos2[i] = tmp_max;
		}	
		DoQuicksort(start_pos1, n1, sort_perm1);   // sort the first positions vec by start
		for(i=0; i<n1; i++)
			tmp_pos[i] = end_pos1[sort_perm1[i]];
		for(i=0; i<n1; i++)
			end_pos1[i] = tmp_pos[i]; 
		DoQuicksort(start_pos2, n2, sort_perm2); // sort the second positions vec by start
		for(i=0; i<n2; i++)
			tmp_pos[i] = end_pos2[sort_perm2[i]];
		for(i=0; i<n2; i++)
			end_pos2[i] = tmp_pos[i]; 
	}
	// print sorted values
	/*
	for(i=0; i<n1; i++)
		printf("sorted edges1: (%lf, %lf)\n", start_pos1[i], end_pos1[i]); 
	for(i=0; i<n2; i++)
		printf("sorted edges2: (%lf, %lf)\n", start_pos2[i], end_pos2[i]); 
*/

//	for(i=0; i<n1; i++)
//		printf("sort_perm1_now %ld\n", sort_perm1[i]);
	

	long prev_ctr2 = 0; // counters following both regions
	double inter_start, inter_end; 
	inter_ctr[0] = 0;

	for(ctr1=0; ctr1<n1; ctr1++) // loop on first interals start positions
	{
		ctr2=-1;
		for(j=prev_ctr2; j<n2; j++)
			if(start_pos1[ctr1] <= end_pos2[j])
			{
				ctr2=j-prev_ctr2+1;
				break;
			}
//		printf("ctr2_is %ld\n", ctr2);
		if(ctr2 > -1) // found a match
		{
			inter_start = MAX(start_pos1[ctr1], start_pos2[prev_ctr2+ctr2-1]);
			inter_end = MIN(end_pos1[ctr1], end_pos2[prev_ctr2+ctr2-1]);
			/**
            printf("before while: start1[%ld] %lf, start2[%ld] %lf\n", ctr1, start_pos1[ctr1], prev_ctr2+ctr2-1, start_pos2[prev_ctr2+ctr2-1]); 
            printf("before while: end1[%ld] %lf, end2[%ld]%lf\n", ctr1, end_pos1[ctr1], prev_ctr2+ctr2-1, end_pos2[prev_ctr2+ctr2-1]); 
			printf("before while inter_start %lf before while inter_end %lf\n", inter_start, inter_end); 
			/**/
			while(inter_start <= inter_end)
			{
				intersect_start_pos[inter_ctr[0]] = inter_start;
				intersect_end_pos[inter_ctr[0]] = inter_end;
				intersect_inds1[inter_ctr[0]] = ctr1;   // New: get also the indexes
				intersect_inds2[inter_ctr[0]] = prev_ctr2+ctr2-1;
				inter_ctr[0] = inter_ctr[0]+1;
//				printf("ctr1 %ld inter_ctr %ld\n", ctr1, inter_ctr[0]); 
				if(prev_ctr2+ctr2 < n2) // make sure we didn't reach the end of second array
				{
					inter_start = MAX(start_pos1[ctr1], start_pos2[prev_ctr2+ctr2]);
					inter_end = MIN(end_pos1[ctr1], end_pos2[prev_ctr2+ctr2]);
				}
				else
				{
					inter_start = 1; inter_end = -1;
				}
//				printf("inside while inter_start %lf inside while inter_end %lf\n", inter_start, inter_end); 
				ctr2=ctr2+1;
			}
//			printf("inter_start2 %lf inter_end2 %lf\n", inter_start, inter_end); 
            prev_ctr2 = MAX(prev_ctr2, prev_ctr2+ctr2-2); // shift previous counter to the right
		}
	}
//	printf("Intersect Length %ld\n", inter_ctr[0]);
//	for(i=0;i<2;i++)
//		printf("start %lf end %lf ind1 %lf ind2 %lf\n", intersect_start_pos[i], intersect_end_pos[i], intersect_inds1[i], intersect_inds2[i]);
	
	/* sorting part: */
	if(!is_sorted) // avoid sort when not needed 
	{
		for(i=0; i<(inter_ctr[0]); i++)
			tmp_pos[i] = sort_perm1[long(intersect_inds1[i])]; 
		for(i=0; i<(inter_ctr[0]); i++)
			intersect_inds1[i] = tmp_pos[i]+1;
		for(i=0; i<(inter_ctr[0]); i++)
			tmp_pos[i] = sort_perm2[long(intersect_inds2[i])]; 
		for(i=0; i<(inter_ctr[0]); i++)
			intersect_inds2[i] = tmp_pos[i]+1;
	}	
	else
	{
		for(i=0; i<(inter_ctr[0]); i++)
			intersect_inds1[i]++;
		for(i=0; i<(inter_ctr[0]); i++)
			intersect_inds2[i]++;
	}
	/**/

/*
	for(i=0; i<n1; i++)
		printf("sort_perm1 %ld\n", sort_perm1[i]);
	for(i=0; i<n2; i++)
		printf("sort_perm2 %ld\n", sort_perm2[i]);
		
	for(i=0; i<inter_ctr[0]; i++)
		printf("inter_ind1 %ld\n", long(intersect_inds1[i]));
	for(i=0; i<inter_ctr[0]; i++)
		printf("inter_ind2 %ld\n", long(intersect_inds2[i]));
	printf("Intersect Length Again %ld\n", inter_ctr[0]);
*/
	if(!is_sorted)
	{
		mxFree(sort_perm1);   // clear memory
		mxFree(sort_perm2);
		mxFree(tmp_pos); 
	}
	return 0;
}


