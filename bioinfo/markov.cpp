#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "general.h"
#include "dna_utils_new.h"
#include "markov.h"







/// Here our problem is the following : given a set of weights, 
/// we want to transfer them into natural numbers .. 
/// We do it by multiplying by a big number and rounding
long transfer_weights(double w[4][4][MAX_L], long m[4][4][MAX_L], double res, long l)
{

	long i, j, k;


	long N = (long)(1.0/res) + 1;


	// multiply everything by N
	for(i = 0; i < 4; i++)
		for(j = 0; j < 4; j++)
			for(k = 0; k < l; k++)
				m[i][j][k] = long ( ceil(N * w[i][j][k]) );


	return N;
}



// return m to the power of d using square and multiply
// This is in order not to trust C's powers.
double my_power(double m, long d)
{

	if(d == 0)
		return 1.0;
	else
	{
		double temp = my_power(m, d >> 1); 
		
		temp = temp*temp;

		if(d&1) // an even number ...
			temp = temp * m;

		return temp;
	}


	
}




// Note : here we also normalize for the promoter length (the trials variable). 
// We assume independent experiments.
// The threshold probability is computed using the generating functions method
double calc_thresh_prob_from_dist(long thresh,  double *dist, long strands, long max_val, long l, long trials)
{
	double thresh_prob = 0.0;

	long i, j;


	if(thresh < max_val)
		for(i = thresh; i < max_val; i++)
			thresh_prob += dist[i];


	// Now normalize for the length ...
	double neg_prob = 1.0-strands*thresh_prob; // if both strands are present we multiply the probability by 2


	double last[MAX_L];
	double temp;
	
	// initilize
	for(i = 0; i < l; i++)
		last[i] = strands*thresh_prob;

	double sum = strands*thresh_prob*l;

	// advance
	for(i = l; i < trials-l+1; i++)  
	{
		// take the new one
		temp = last[0] - strands*thresh_prob*last[l-1];

		// "shift" the register 
		for(j = l-1; j >= 1; j--)
			last[j] = last[j-1];
		last[0] = temp;
		sum += temp;
	}





//	return 1.0 - my_power(neg_prob, trials);
	
	return sum; //last[0];
//	return thresh_prob;

}




// Calculate the probability to be above a certain threshold
// What is the difference from the previous function ???
double calc_thresh_prob(long m[4][4][MAX_L], long thresh, long l, double *probs, double *dist, long *max_val, long *min_val)
{
	long i, j, k, r;

	double markov_probs[1UL << 13];

	long max_score = 0;
	long min_score = 0;
	long range;
	
	long prev[16];
	long curr[16];

	double *v_prev[16];
	double *v_curr[16];
	
	long rm[4][4][MAX_L];
	long rthresh;

	double total_prob = 0;


	FILE *dist_f;          // distribution file
	
	char *dist_p = "score_dist.txt";			// 	distribution file name



	probs_to_markov_probs(probs, markov_probs, 2); // currently only dim 2 is supported ..


	for(i = 0; i < 4; i++)  
		for(j = 0; j < 4; j++)
			for(k = 0; k < l; k++)  // advance the layers ..
				if(min_score > m[i][j][k])
					min_score = m[i][j][k];



	// generate the mmm's 
	for(i = 0; i < 4; i++)  
		for(j = 0; j < 4; j++)
			for(k = 0; k < l; k++)  // advance the layers ..
				rm[i][j][k] = m[i][j][k] - min_score;
	rthresh = thresh - l*min_score;

	
	*min_val = l*min_score;


	// First find the maximal N possible, using viterbi ...	
	for(i = 0; i < 16; i++)
		prev[i] = 0;

	for(k = 0; k < l; k++)  // advance the layers ...
	{
		for(i = 0; i < 16; i++)  // go over the 16 "new nodes"
		{
			curr[i] = prev[i >> 2]; // + m[i >> 2][i&0x3][k];  // try the first ..

			for(j = 1; j < 4; j++)
				if(prev[(j << 2) + (i >> 2)] > curr[i])
					curr[i] = prev[(j << 2) + (i >> 2)];

			curr[i] += rm[i >> 2][i&0x3][k];

		}

		// change curr to be prev
		for(i = 0; i < 16; i++)
			prev[i] = curr[i];

	}


	// Take the maximum ...
	for(i = 0; i < 16; i++)
		if(max_score < curr[i])
			max_score = curr[i];





	// Second find the minimal N possible, using viterbi ...	
	for(i = 0; i < 16; i++)
		prev[i] = 0;

	for(k = 0; k < l; k++)  // advance the layers ...
	{
		for(i = 0; i < 16; i++)  // go over the 16 "new nodes"
		{
			curr[i] = prev[i >> 2]; // + m[i >> 2][i&0x3][k];  // try the first ..

			for(j = 1; j < 4; j++)
				if(prev[(j << 2) + (i >> 2)] < curr[i])
					curr[i] = prev[(j << 2) + (i >> 2)];

			curr[i] += rm[i >> 2][i&0x3][k];

		}

		// change curr to be prev
		for(i = 0; i < 16; i++)
			prev[i] = curr[i];

	}


	// Take the minimum ...
	for(i = 0; i < 16; i++)
		if(min_score > curr[i])
			min_score = curr[i];

//	printf("The Maximum and Minimum are : %ld   %ld\n", max_score, min_score);


/**	// changing it to include l ..
	min_score = min_score * -1;
	min_score = ((min_score/l)+1)*l;

//	printf("NOW MIN IS %ld\n");


	// generate the mmm's 
	for(i = 0; i < 4; i++)  
		for(j = 0; j < 4; j++)
			for(k = 0; k < l; k++)  // advance the layers ..
				rm[i][j][k] = m[i][j][k] + min_score/l;
	rthresh = thresh + min_score;
***/


/**	
	printf("The long meights ..\n");
	printf("  AA   AT   AC   AG   TA   TT   TC   TG   CA   CT   CC   CG   GA   GT   GC   GG\n");
	// print the long weights 
	for(k = 0; k < l; k++)
	{
		for(i = 0; i < 4; i++)
			for(j = 0; j < 4; j++)
				printf("%4ld ", m[i][j][k]);
	
				
		printf("\n");

	}

	printf("The long rrrrrmeights ..\n");
	printf("  AA   AT   AC   AG   TA   TT   TC   TG   CA   CT   CC   CG   GA   GT   GC   GG\n");
	// print the long weights 
	for(k = 0; k < l; k++)
	{
		for(i = 0; i < 4; i++)
			for(j = 0; j < 4; j++)
				printf("%4ld ", rm[i][j][k]);
	
				
		printf("\n");

	}
**/




	range = max_score/*+min_score*/+1;


	// Now the score will always be between 0 and range ...	


	// We need max - min
	for(i = 0; i < 16; i++)
	{
		v_prev[i] = new double[range];
		v_curr[i] = new double[range];

		for(r = 0; r < range; r++)
		{
			v_prev[i][r] = 0.0;
			v_curr[i][r] = 0.0;
		}
		
		v_prev[i][0] = 1.0/16.0;  // start with zero in the score

	}	



	// Now do the serious staff ....
	// here k=0 , we use probs ..
	for(i = 0; i < 16; i++)  // go over the 16 "new nodes"
	{
			
		// go over the nodes that can lead to our node
		for(j = 0; j < 4; j++)
			for(r = 0; r < range; r++)
				if(v_prev[(j << 2) + (i >> 2)][r] != 0)  // we have a positive probability
				{
					//printf("new cell %ld\n", r + rm[i >> 2][i&0x3][k]);
					v_curr[i][ r + rm[i >> 2][i&0x3][0] ] += 
						v_prev[(j << 2) + (i >> 2)][r] * probs[i]; 
				}

	}



	// set the current to zero again ..
	for(i = 0; i < 16; i++)
		for(r = 0; r < range; r++)	
		{
			v_prev[i][r] = v_curr[i][r];
			v_curr[i][r] = 0.0;
		}




	// Now do the serious staff .... k> 0
	for(k = 1; k < l; k++)  // advance the layers ...
	{
		for(i = 0; i < 16; i++)  // go over the 16 "new nodes"
		{
			
			// go over the nodes that can lead to our node
			for(j = 0; j < 4; j++)
				for(r = 0; r < range; r++)
					if(v_prev[(j << 2) + (i >> 2)][r] != 0)  // we have a positive probability
					{
						//printf("new cell %ld\n", r + rm[i >> 2][i&0x3][k]);
						v_curr[i][ r + rm[i >> 2][i&0x3][k] ] += 
							v_prev[(j << 2) + (i >> 2)][r] * markov_probs[i]; // >> 2][i&0x3]/*[k]*/;
					}

		}

/**
		// change curr to be prev
		for(i = 0; i < 16; i++)
		{
		//	prev[i] = curr[i];
			SWAP(v_prev[i], v_curr[i], v_temp);
		}
**/


		// set the current to zero again ..
		for(i = 0; i < 16; i++)
			for(r = 0; r < range; r++)	
			{
				v_prev[i][r] = v_curr[i][r];
				v_curr[i][r] = 0.0;
			}
	}





	

	// Now collect the probabilities 
	total_prob = 0;
	for(i = 0; i < 16; i++)
		for(r = 0; r < range; r++)
			total_prob += v_prev[i][r];

	
	total_prob = 0;
	for(i = 0; i < 16; i++)
		for(r = MAX(rthresh, 0); r < range; r++)
		{
			total_prob += v_prev[i][r];

		}


	
	*max_val = range;
	for(r = 0; r < range; r++)
		dist[r] = v_prev[0][r];

	// use dist to accumelate ...
	for(i = 1; i < 16; i++)
		for(r = 0; r < range; r++)
			dist[r] /*v_prev[0][r]*/ += v_prev[i][r];



	// open the file to enable writing to it
	dist_f = fopen(dist_p, "w+");

	// get back to the m's ...
	for(r = 0; r < range; r++)
		fprintf(dist_f, "%ld %lf\n", r + l*min_score, dist[r]); 


	fclose(dist_f);


	// Start by the first 

/**/


	// We need max - min
	for(i = 15; i >= 0; i--)
	{
		delete 	v_prev[i];
		delete 	v_curr[i];
	}	
/**/


	return total_prob;


}




