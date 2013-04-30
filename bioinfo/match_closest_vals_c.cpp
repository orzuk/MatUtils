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
/// Find closest values in a to b - 07.10.2009
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
				  int nrhs, const mxArray *prhs[])
{
	long i, j, a_len, b_len;	
	double *a, *b;
	double *in, *inn, *out;

	/* Check for proper number of arguments. Should be Two 2 */
	if(nrhs != 2) {
		mexErrMsgTxt("Usage: [dat] [K]");
	}  
	in = mxGetPr(prhs[0]);   // Get the first values vector matrix
	a_len = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));    // Get the number of sequences (genes)
	a = new double[a_len];
	for(i = 0; i < a_len; i++)
		a[i] = in[i];  // copy a values
	in = mxGetPr(prhs[1]);   // Get the first values vector matrix
	b_len = MAX(mxGetM(prhs[1]), mxGetN(prhs[1]));    // Get the number of sequences (genes)
	b = new double[b_len];
	for(i = 0; i < b_len; i++)
		b[i] = in[i];  // copy b values


	double *closest_ind = new double[a_len];
	double *closest_dist = new double[a_len];
	match_closest_vals(closest_ind, closest_dist, a, b, a_len, b_len); // call function
	
	// The two outputs of the function 
	plhs[0] = mxCreateDoubleMatrix(a_len, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(a_len, 1, mxREAL);

	out = mxGetPr(plhs[0]); // output the inds array 
	for(i=0; i<a_len; i++)
		out[i] = closest_ind[i];
	out = mxGetPr(plhs[1]); // output the distances array 
	for(i=0; i<a_len; i++)
		out[i] = closest_dist[i];

	delete a, b, closest_ind, closest_dist; 	// free memory

}
//#endif // MEX_TO_MATLAB


// Find for each element of a the closest element of b (provide speedup for matlab function)
// We assume that they are already sorted and allocated 

/*** Matlab code here: 
if(isempty(b))
    closest_ind = []; closest_dist = []; return; % nothing to compare to 
end
[a_sorted sort_perm_a] = sort(a); n = length(a);
[b_sorted sort_perm_b] = sort(b); m = length(b);
inv_sort_perm_a = inv_perm(sort_perm_a); 
% inv_sort_perm_b = inv_perm(sort_perm_b); 
closest_ind = zeros(n,1); closest_dist = zeros(n,1); 
ctr = 1;
for i=1:n
    f = find(a_sorted(i) < b_sorted(ctr:end), 1);
    if(isempty(f)) % this means that this value of a is bigger than all values of b     
        closest_ind(i:n) = m;
        closest_dist(i:n) = a_sorted(i:n) - b_sorted(m); 
        break; % get out of loop
    else
        if(f == 1)
            closest_ind(i) = ctr;
        else
            [dummy closest_ind(i)] = min(abs(b_sorted(ctr+f-2:ctr+f-1) - a_sorted(i)));
            closest_ind(i) =  closest_ind(i) + f-2 + ctr-1;
        end
        closest_dist(i) = a_sorted(i) - b_sorted(closest_ind(i)); 
        ctr = closest_ind(i); 
    end
end
% apply inverse permutations to return to original
closest_ind = sort_perm_b(closest_ind(inv_sort_perm_a)); 
closest_dist = closest_dist(inv_sort_perm_a); 
***/
long match_closest_vals(double *closest_ind, double *closest_dist, double *a, double *b, long a_len, long b_len)
{
	long i,j,f;
	long ctr = 0;

	for(i=0; i<a_len; i++)
	{
		f=-1;
		for(j=ctr; j<b_len; j++)
		{
			if(a[i] < b[j])
			{
				f=j-ctr;
				break;
			}
		}
		if(f==-1) // couldn't find and element
		{
			for (j=i; j<a_len; j++)
			{
				closest_ind[j] = b_len;
				closest_dist[j] = a[j] - b[b_len-1]; 
			}
		    break; // get out of loop
		}
		else
		{
			if(f==0)
				closest_ind[i] = ctr;
			else
			{
				if(abs(b[ctr+f-1] - a[i]) <= abs(b[ctr+f] - a[i]))
					closest_ind[i]=0;
				else
					closest_ind[i]=1;
				closest_ind[i] =  closest_ind[i] + f-1 + ctr; // compensate two - one for closest_ind and one for f
			}
			closest_dist[i] = a[i] - b[long(closest_ind[i])];
			ctr = closest_ind[i];
		}
		closest_ind[i] = closest_ind[i]+1; // convert C 0-start to Matlab 1-start
	}
	return 0;
}

