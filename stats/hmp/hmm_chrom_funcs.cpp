#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "general.h"
#include "hmm_chrom_funcs.h"

#undef DEBUG


/////////////////////////////////////////////////////////////////////
// Matrix multiplication routine 
// A is m*n
// B is k*m
// Output C = A*B is r*n. We assume C is allocated 
// we assume the maximal dimentions are MAX_X_VALS 
//
//
long MatrixMultiply(double A[MAX_X_VALS][MAX_X_VALS], long m, long n, double B[MAX_X_VALS][MAX_X_VALS], long r, 
					double C[MAX_X_VALS][MAX_X_VALS])
{
	long i, j, k;

	
	// First init C to zero
	for(i = 0; i < r; i++)
		for(j = 0; j < n; j++)
			C[i][j] = 0; 


	// Now compute C 
	for( i = 0; i < n; i++ )   // Go over the rows of A
         for( j = 0; j < r; j++ )
             for( k = 0; k < m; k++ )   // go over elements of the line
				C[i][j] += A[i][k] * B[k][j];


	return 0; 


}
/////////////////////////////////////////////////////////////////////
long MatrixTranspose(double A[MAX_X_VALS][MAX_X_VALS], long m, long n, 
					 double A_t[MAX_X_VALS][MAX_X_VALS])

{
	long i, j;

	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			A_t[j][i] = A[i][j];
	return 0;
}




// Get the stationary distribution for a markov (stochastic) matrix A
long MatrixStationaryVec(double A[MAX_X_VALS][MAX_X_VALS], long n, double pi[MAX_X_VALS])
{
	long i,j;
	double sum_pi=0;

	double A_t[MAX_X_VALS][MAX_X_VALS];

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			A_t[i][j] = A[i][j];

	for(i=0;i<n;i++)
		A_t[i][i] -= 1; // subtract the identity matrix

	for(i=0;i<n;i++)
		for(j=0;j<n-1;j++)
			A_t[j][i] -= A_t[n-1][i]; // subtract the last column

	for(i=0;i<n;i++)
		A_t[n-1][i] *= (-1); // transfer sides to the B vector

	GSolve(A_t, n-1, pi);

	for(i=0;i<n-1;i++)
		sum_pi += pi[i];
	pi[n-1] = 1-sum_pi;

	
	return 0;
}



/*
   Solve a system of n equations in n unknowns using Gaussian Elimination
   Solve an equation in matrix form Ax = b
   The 2D array a is the matrix A with an additional column b.
   This is often written (A:b)

   A0,0    A1,0    A2,0    ....  An-1,0     b0
   A0,1    A1,1    A2,1    ....  An-1,1     b1
   A0,2    A1,2    A2,2    ....  An-1,2     b2
   :       :       :             :          :
   :       :       :             :          :
   A0,n-1  A1,n-1  A2,n-1  ....  An-1,n-1   bn-1

   The result is returned in x, otherwise the function returns FALSE
   if the system of equations is singular.
*/
long GSolve(double a[MAX_X_VALS][MAX_X_VALS],long n,double x[MAX_X_VALS])
{
   long i,j,k,maxrow;
   double tmp;

   for (i=0;i<n;i++) 
   {
      /* Find the row with the largest first value */
      maxrow = i;
      for (j=i+1;j<n;j++) {
         if (ABS(a[i][j]) > ABS(a[i][maxrow]))
            maxrow = j;
      }

      /* Swap the maxrow and ith row */
      for (k=i;k<n+1;k++) {
         tmp = a[k][i];
         a[k][i] = a[k][maxrow];
         a[k][maxrow] = tmp;
      }

      /* Singular matrix? */
      if (ABS(a[i][i]) < EPSILON)
         return(FALSE);

      /* Eliminate the i-th element of the j-th row */
      for (j=i+1;j<n;j++) {
         for (k=n;k>=i;k--) {
            a[k][j] -= a[k][i] * a[i][j] / a[i][i];
         }
      }
   }

   /* Do the back substitution */
   for (j=n-1;j>=0;j--) {
      tmp = 0;
      for (k=j+1;k<n;k++)
         tmp += a[k][j] * x[k];
      x[j] = (a[n][j] - tmp) / a[j][j];
   }

   return(TRUE);
}










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

// inner split routine: overloaded for float
long split(float *vals, long first, long last, long *indexes)
{
    float    x, temp;
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


// recursive quicksort routine: overloaded for float
long quicksort(float *vals, long first, long last, long *indexes)
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
		indexes[i] = i;  // init: indexes are i.d. permutation 		
	quicksort(vals, 0, len-1, indexes);
	return 0; 
}

// outer  call for quicksort: overloaded function for floats 
long DoQuicksort(float *vals, long len, long *indexes)
{
	long i;
	for(i = 0; i < len; i++)
		indexes[i] = i;  // init: indexes are i.d. permutation 		
	quicksort(vals, 0, len-1, indexes);
	return 0; 
}


// order an area by a permutation
long DoOrder(double *vals, long len, long *indexes)
{
	long i;
	double *tmp_arr = new double[len];
	for(i = 0; i < len; i++)
		tmp_arr[i] = vals[indexes[i]];
	for(i = 0; i < len; i++)
		vals[i] = tmp_arr[i];
	delete tmp_arr;
	return 0; 
}

// order an area by a permutation: overloaded for float
long DoOrder(float *vals, long len, long *indexes)
{
	long i;
	float *tmp_arr = new float[len];
	for(i = 0; i < len; i++)
		tmp_arr[i] = vals[indexes[i]];
	for(i = 0; i < len; i++)
		vals[i] = tmp_arr[i];
	delete tmp_arr;
	return 0; 
}



// Print a vector of integers
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

// Print a vector of doubles
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

// generate a random double in the interval [0,1]
double randr()
{
	return double(rand())/double((RAND_MAX+1.0));
}


// Generate a geometric variable with parameter p. p is 'success' probability
long Geometric(double p)
{
	long failures = 0; 

	// we assume we don't have too large gaps so actually we get a truncated geometric distribution
	while((randr() > p) && (failures < MAX_POWER_OF_TRANS_MATRIX-1))   
		failures++;

	return failures;

}



// Generates a gaussian random variable with mean mu and standard deviation sigma.
// Method used is the Box-Miller transform 
double Gaussian(double mu,double sigma)

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
		return(mu+v1*r*sigma);
   }
   else 
   {
     x = t;
     t = 0.0;
     return(mu+x*sigma);
   }
 }




// Here we give some Gaussian Mixture routines: learning, classification, parameter estimation etc.
// All these are implemented in Matlab but we want to get high-performance so we put them here.
// We implement seperately one dimensional and multi-dimensional models.
double MixtureOfGaussiansGivenInit(double *x, long num_points, long num_of_Gaussians, long num_of_iterations, 
		 double *init_p, double *init_m, double *init_s, 
		 double *p, double *m, double *s, double *log_like)
{
	long i, j, k, g, itt, iter;
	double KL, tmp_KL;

	if(num_of_Gaussians == 1) // Deal first with the trivial case of one gaussian
	{
		m[0] = 0;
		for(i=0; i<num_points; i++)
		{
			m[0] += x[i];
			s[0] += (x[i]*x[i]);
		}
		m[0] = m[0] / num_points;
		s[0] = s[0] / num_points - m[0]*m[0];
		p[0] = 1;
	    (*log_like) = -num_points*log(SQRT_2PI*s[0]);
		for(i=0; i<num_points; i++)
			(*log_like) -= ( (x[i]-m[0])*(x[i]-m[0]) / (2*s[0]*s[0]) );
	}
	else
	{
	    double TOL=0.000000001; // tolerance for score improvement
		double maxKL=-100000000000.0; // start with very very low log-likelihood
	    double OldKL = -9999999999999999.9;
//		double p_prior_mult = 0; 

		double sigma[MAX_NUM_GAUSSIANS];
		double mu[MAX_NUM_GAUSSIANS];
		double prior[MAX_NUM_GAUSSIANS];
		double prev_prior[MAX_NUM_GAUSSIANS];
		
		double *pr[MAX_NUM_GAUSSIANS];
		double *z[MAX_NUM_GAUSSIANS];
		double dummy;

		double sum_z[MAX_NUM_GAUSSIANS];
		long m_inds[MAX_NUM_GAUSSIANS];

		double *p_prior_mult = new double[num_points];

		for(g=0; g<num_of_Gaussians; g++)
		{
			pr[g] = new double[num_points]; 
			z[g] = new double[num_points]; 
//			y[g] = new double[num_points]; 
			sigma[g] = init_s[g];
			mu[g] = init_m[g];
			prior[g] = init_p[g];
		}

		
		for (itt=0; itt<num_of_iterations; itt++) // Start EM iterations - this is the heavy loop
		{
			/*******************  E-step: update counts ********************************************/
	        for(g=0; g<num_of_Gaussians; g++) 			
			{
				if(sigma[g] > 0)
					for(i=0; i<num_points; i++)
						pr[g][i] = exp(- ((x[i]-mu[g])*(x[i]-mu[g]) / (sigma[g]*sigma[g]*2))  ) / (SQRT_2PI*sigma[g]);
				else
					for(i=0; i<num_points; i++)
			            pr[g][i]=1;
			}
			for(i=0; i<num_points; i++)
				p_prior_mult[i] = pr[0][i] * prior[0];
			for (g=1; g<num_of_Gaussians; g++)
				for(i=0; i<num_points; i++)
					p_prior_mult[i] += pr[g][i] * prior[g];
			for (g=0; g<num_of_Gaussians; g++)
				for(i=0; i<num_points; i++)
					z[g][i] = pr[g][i] * prior[g] / p_prior_mult[i];
			/*******************  Ended E-step *********************************************************/
			/*******************  M-step: update parameters ********************************************/
			for(g=0; g<num_of_Gaussians; g++) 	         
			{
				sum_z[g]=0;
				for(i=0; i<num_points; i++)
					sum_z[g] += z[g][i];
				if(sum_z[g] == 0)
					sum_z[g] = 1;
				mu[g]=0;    // update mu
				for(i=0; i<num_points; i++)
					mu[g] += z[g][i] * x[i];
				mu[g] /= sum_z[g];
				sigma[g]=0; 				// update sigma AFTER mu - with the new mu
				for(i=0; i<num_points; i++)
					sigma[g] += (z[g][i] * x[i] * x[i]);
				sigma[g] = sqrt(sigma[g] / sum_z[g] - mu[g]*mu[g]);
//				for(i=0; i<num_points; i++)
//					sigma[g] += z[g][i] * (x[i]-mu[g]) * (x[i]-mu[g]);
//				sigma[g] = sqrt(sigma[g] / sum_z[g]);
				sigma[g] = MAX(sigma[g], TOL); //  Regularization: Do not allow too small sigmas
				prev_prior[g] = prior[g]; prior[g] = sum_z[g] / num_points;
			/*******************  Ended M-step *********************************************************/
			}



	        // Calculate the score, and check if improvement is not negligible
			KL=0.0;
			for(i=0; i<num_points; i++)
			{
//				tmp_KL = prior[0] / (SQRT_2PI*sigma[0]) * exp ( -(x[i]-mu[0])*(x[i]-mu[0]) / (2*sigma[0]*sigma[0]) );
//				for(g=1; g<num_of_Gaussians; g++)
//					tmp_KL += ( prior[g] / (SQRT_2PI*sigma[g]) * exp ( -(x[i]-mu[g])*(x[i]-mu[g]) / (2*sigma[g]*sigma[g]) ) );
													  

				// get the likelihood of previous round - a bit faster (is this incorrect??) 
				tmp_KL = prev_prior[0] * pr[0][i];
				for(g=1; g<num_of_Gaussians; g++)
					tmp_KL += prev_prior[g] * pr[g][i];				
				KL += log(tmp_KL); // log(tmp_KL);
			}
	        if(KL - OldKL < TOL)
			{
		        iter = itt;
			    break;
			}
			OldKL = KL;
		}  // end looping on number of iterations



		// Check for improvement in score (with respect to previous starting
		// points) and update parameters
/***/	if(KL>maxKL) /***/
		{
		    (*log_like) = KL; maxKL=KL;
			for(g=0; g<num_of_Gaussians; g++)
			{
				p[g]=prior[g]; s[g]=sigma[g]; m[g]=mu[g];
			}
	    }	
		

	    // Sort the means in increasing order		
		long qsort_res = DoQuicksort(m, num_of_Gaussians, m_inds);
		DoOrder(p, num_of_Gaussians, m_inds);
		DoOrder(s, num_of_Gaussians, m_inds); 

		(*log_like) = maxKL; 	    // return also the best log-likelihood

		delete p_prior_mult; 
		for(g=0; g<num_of_Gaussians; g++)
		{
			delete pr[g];
			delete z[g];
//			delete y[g];
		}
    }
	return double(iter); 

///////////////////////////////////////////////////////////////
/// Start Matlab Code ...
///////////////////////////////////////////////////////////////
/***
	N=length(x); % # data points
if(num_of_Gaussians == 1) % Deal first with the trivial case of one gaussian
    M = mean(x); S = std(x); P = 1;
    LogLike = -N*log(sqrt(2*pi)*S) - sum((x-M).^2 ./ (2*S^2));
else

    p=ones(N,num_of_Gaussians)/num_of_Gaussians; % start with uniform probs.
    z=p; % just allocate memory to save some time .. 
    
    two_pi_sqrt = sqrt(2*pi); % A constant to save time
    TOL=0.000000001; % tolerance for score improvement
    maxKL=-1000000000; % start with very very low log-likelihood

    OldKL = -9999999999999999;
    sigma=INIT_S; miu = INIT_M; prior = INIT_P;
    for itt=1:num_of_iterations % Start EM iterations - this is the heavy loop
        % E-step: update counts
        for m=1:num_of_Gaussians
            if(sigma(m))
                p(:,m)=exp(- (((x-miu(m))./(sigma(m))).^2) ./ 2 ) ./ (two_pi_sqrt.*sigma(m));
            else
                p(:,m)=1;
            end
        end
        for m=1:num_of_Gaussians
            z(:,m)=(p(:,m).*prior(m))./((p*prior'));
        end

        % M-step: update parameters
        sum_z=sum(z); sum_z(sum_z==0)=1;
        sigma=sqrt(sum(z.*(repmat(x,1,num_of_Gaussians)-repmat(miu,N,1)).^2)./sum_z);          %x is a column vector
        sigma = max(sigma,TOL); % Regularization: Do not allow too small sigmas
        miu=(z'*x)'./sum_z;
        prior=sum_z/N;

        % Calculate the score, and check if improvement is not negligible
        for m=1:num_of_Gaussians
            y(m,:)=prior(m)*1/(two_pi_sqrt*sigma(m))*exp(-( (x-miu(m))./sigma(m)).^2./2);
        end
        KL=sum(log(sum(y,1)));
        if(KL - OldKL < TOL)
            iter = itt
            break;
        end
        OldKL = KL;
    end

    % Check for improvement in score (with respect to previous starting
    % points) and update parameters
    if(KL>maxKL)
        LogLike = KL % print current score to the screen
        maxKL=KL;
        P=prior; S=sigma; M=miu;
    end


    % Sort the means in increasing order
    [M IndM] = sort(M);
    P = P(IndM); S = S(IndM);
    % return also the best log-likelihood
    LogLike = maxKL;

end
***/
  ///////////////////////////////////////////////////////////////
/// Ended Matlab Code ...
///////////////////////////////////////////////////////////////

}



// Here give the same MoG function, but for single rather than double !!! (overload) 
// Here we give some Gaussian Mixture routines: learning, classification, parameter estimation etc.
// All these are implemented in Matlab but we want to get high-performance so we put them here.
// We implement seperately one dimensional and multi-dimensional models.
float MixtureOfGaussiansGivenInitSingle(float *x, long num_points, long num_of_Gaussians, long num_of_iterations, 
		 float *init_p, float *init_m, float *init_s, 
		 float *p, float *m, float *s, float *log_like)
{
	long i, j, k, g, itt, iter;
	float KL, tmp_KL;

	if(num_of_Gaussians == 1) // Deal first with the trivial case of one gaussian
	{
		m[0] = 0;
		for(i=0; i<num_points; i++)
		{
			m[0] += x[i];
			s[0] += (x[i]*x[i]);
		}
		m[0] = m[0] / num_points;
		s[0] = s[0] / num_points - m[0]*m[0];
		p[0] = 1;
	    (*log_like) = -num_points*log(SQRT_2PI*s[0]);
		for(i=0; i<num_points; i++)
			(*log_like) -= ( (x[i]-m[0])*(x[i]-m[0]) / (2*s[0]*s[0]) );
	}
	else
	{
	    float TOL=0.000000001; // tolerance for score improvement
		float maxKL=-100000000000.0; // start with very very low log-likelihood
	    float OldKL = -9999999999999999.9;
//		double p_prior_mult = 0; 

		float sigma[MAX_NUM_GAUSSIANS];
		float mu[MAX_NUM_GAUSSIANS];
		float prior[MAX_NUM_GAUSSIANS];
		float prev_prior[MAX_NUM_GAUSSIANS];
		
		float *pr[MAX_NUM_GAUSSIANS];
		float *z[MAX_NUM_GAUSSIANS];
		float dummy;

		float sum_z[MAX_NUM_GAUSSIANS];
		long m_inds[MAX_NUM_GAUSSIANS];

		float *p_prior_mult = new float[num_points];

		for(g=0; g<num_of_Gaussians; g++)
		{
			pr[g] = new float[num_points]; 
			z[g] = new float[num_points]; 
//			y[g] = new float[num_points]; 
			sigma[g] = init_s[g];
			mu[g] = init_m[g];
			prior[g] = init_p[g];
		}

		
		for (itt=0; itt<num_of_iterations; itt++) // Start EM iterations - this is the heavy loop
		{
			/*******************  E-step: update counts ********************************************/
	        for(g=0; g<num_of_Gaussians; g++) 			
			{
				if(sigma[g] > 0)
					for(i=0; i<num_points; i++)
						pr[g][i] = exp(- ((x[i]-mu[g])*(x[i]-mu[g]) / (sigma[g]*sigma[g]*2))  ) / (SQRT_2PI*sigma[g]);
				else
					for(i=0; i<num_points; i++)
			            pr[g][i]=1;
			}
			for(i=0; i<num_points; i++)
				p_prior_mult[i] = pr[0][i] * prior[0];
			for (g=1; g<num_of_Gaussians; g++)
				for(i=0; i<num_points; i++)
					p_prior_mult[i] += pr[g][i] * prior[g];
			for (g=0; g<num_of_Gaussians; g++)
				for(i=0; i<num_points; i++)
					z[g][i] = pr[g][i] * prior[g] / p_prior_mult[i];
			/*******************  Ended E-step *********************************************************/
			/*******************  M-step: update parameters ********************************************/
			for(g=0; g<num_of_Gaussians; g++) 	         
			{
				sum_z[g]=0;
				for(i=0; i<num_points; i++)
					sum_z[g] += z[g][i];
				if(sum_z[g] == 0)
					sum_z[g] = 1;
				mu[g]=0;    // update mu
				for(i=0; i<num_points; i++)
					mu[g] += z[g][i] * x[i];
				mu[g] /= sum_z[g];
				sigma[g]=0; 				// update sigma AFTER mu - with the new mu
				for(i=0; i<num_points; i++)
					sigma[g] += (z[g][i] * x[i] * x[i]);
				sigma[g] = sqrt(sigma[g] / sum_z[g] - mu[g]*mu[g]);
//				for(i=0; i<num_points; i++)
//					sigma[g] += z[g][i] * (x[i]-mu[g]) * (x[i]-mu[g]);
//				sigma[g] = sqrt(sigma[g] / sum_z[g]);
				sigma[g] = MAX(sigma[g], TOL); //  Regularization: Do not allow too small sigmas
				prev_prior[g] = prior[g]; prior[g] = sum_z[g] / num_points;
			/*******************  Ended M-step *********************************************************/
			}



	        // Calculate the score, and check if improvement is not negligible
			KL=0.0;
			for(i=0; i<num_points; i++)
			{
//				tmp_KL = prior[0] / (SQRT_2PI*sigma[0]) * exp ( -(x[i]-mu[0])*(x[i]-mu[0]) / (2*sigma[0]*sigma[0]) );
//				for(g=1; g<num_of_Gaussians; g++)
//					tmp_KL += ( prior[g] / (SQRT_2PI*sigma[g]) * exp ( -(x[i]-mu[g])*(x[i]-mu[g]) / (2*sigma[g]*sigma[g]) ) );
													  

				// get the likelihood of previous round - a bit faster (is this incorrect??) 
				tmp_KL = prev_prior[0] * pr[0][i];
				for(g=1; g<num_of_Gaussians; g++)
					tmp_KL += prev_prior[g] * pr[g][i];				
				KL += log(tmp_KL); // log(tmp_KL);
			}
	        if(KL - OldKL < TOL)
			{
		        iter = itt;
			    break;
			}
			OldKL = KL;
		}  // end looping on number of iterations



		// Check for improvement in score (with respect to previous starting
		// points) and update parameters
/***/	if(KL>maxKL) /***/
		{
		    (*log_like) = KL; maxKL=KL;
			for(g=0; g<num_of_Gaussians; g++)
			{
				p[g]=prior[g]; s[g]=sigma[g]; m[g]=mu[g];
			}
	    }	
		

	    // Sort the means in increasing order		
		long qsort_res = DoQuicksort(m, num_of_Gaussians, m_inds);
		DoOrder(p, num_of_Gaussians, m_inds);
		DoOrder(s, num_of_Gaussians, m_inds); 

		(*log_like) = maxKL; 	    // return also the best log-likelihood

		delete p_prior_mult; 
		for(g=0; g<num_of_Gaussians; g++)
		{
			delete pr[g];
			delete z[g];
//			delete y[g];
		}
    }
	return float(iter); 

///////////////////////////////////////////////////////////////
/// Start Matlab Code ...
///////////////////////////////////////////////////////////////
/***
	N=length(x); % # data points
if(num_of_Gaussians == 1) % Deal first with the trivial case of one gaussian
    M = mean(x); S = std(x); P = 1;
    LogLike = -N*log(sqrt(2*pi)*S) - sum((x-M).^2 ./ (2*S^2));
else

    p=ones(N,num_of_Gaussians)/num_of_Gaussians; % start with uniform probs.
    z=p; % just allocate memory to save some time .. 
    
    two_pi_sqrt = sqrt(2*pi); % A constant to save time
    TOL=0.000000001; % tolerance for score improvement
    maxKL=-1000000000; % start with very very low log-likelihood

    OldKL = -9999999999999999;
    sigma=INIT_S; miu = INIT_M; prior = INIT_P;
    for itt=1:num_of_iterations % Start EM iterations - this is the heavy loop
        % E-step: update counts
        for m=1:num_of_Gaussians
            if(sigma(m))
                p(:,m)=exp(- (((x-miu(m))./(sigma(m))).^2) ./ 2 ) ./ (two_pi_sqrt.*sigma(m));
            else
                p(:,m)=1;
            end
        end
        for m=1:num_of_Gaussians
            z(:,m)=(p(:,m).*prior(m))./((p*prior'));
        end

        % M-step: update parameters
        sum_z=sum(z); sum_z(sum_z==0)=1;
        sigma=sqrt(sum(z.*(repmat(x,1,num_of_Gaussians)-repmat(miu,N,1)).^2)./sum_z);          %x is a column vector
        sigma = max(sigma,TOL); % Regularization: Do not allow too small sigmas
        miu=(z'*x)'./sum_z;
        prior=sum_z/N;

        % Calculate the score, and check if improvement is not negligible
        for m=1:num_of_Gaussians
            y(m,:)=prior(m)*1/(two_pi_sqrt*sigma(m))*exp(-( (x-miu(m))./sigma(m)).^2./2);
        end
        KL=sum(log(sum(y,1)));
        if(KL - OldKL < TOL)
            iter = itt
            break;
        end
        OldKL = KL;
    end

    % Check for improvement in score (with respect to previous starting
    % points) and update parameters
    if(KL>maxKL)
        LogLike = KL % print current score to the screen
        maxKL=KL;
        P=prior; S=sigma; M=miu;
    end


    % Sort the means in increasing order
    [M IndM] = sort(M);
    P = P(IndM); S = S(IndM);
    % return also the best log-likelihood
    LogLike = maxKL;

end
***/
  ///////////////////////////////////////////////////////////////
/// Ended Matlab Code ...
///////////////////////////////////////////////////////////////

}




// A simple function which calculates the joint log-likelihood of the data given the model.
// Nothing tricky here - just the trivial algorithm. We assume for now that Y dim is one
double ComputeJointHMMLogLikelihood(hmm_model *hmm, hmm_data *data)
{
	double loglike_score = 0.0; // This is what we return
	long t;
	double pi = 3.1415927; 

	// First get the score for the X's 
	for(t = 0; t < data->seq_len-1; t++)
	{
		loglike_score += log( hmm->M[data->x_vec[t]][data->x_vec[t+1]] ) - 
			(data->y_vec[t] - hmm->MU[data->x_vec[t]][0]) * (data->y_vec[t] - hmm->MU[data->x_vec[t]][0]) / 
			(2*hmm->SIGMA[data->x_vec[t]][0]*hmm->SIGMA[data->x_vec[t]][0]) - 
			0.5*log(2*pi) - log(hmm->SIGMA[data->x_vec[t]][0]); 
	}

	return loglike_score;
}



// A simple function which calculates the joint log-likelihood of the data given the model.
// Nothing tricky here - just the trivial algorithm. This is for the SNPs model.
// We assume for now that Y dim is one
double ComputeJointHMMLogLikelihoodSNPs(hmm_model *hmm, hmm_data *data,
										double *x_loglike, double *place_x_loglike, double *y_loglike)
{

	double loglike_score = 0.0; // This is what we return
	long t;
	long cur_x_alpha_geno, cur_x_alpha_copy, cur_x_beta_geno, cur_x_beta_copy;
	long next_x_alpha_geno, next_x_alpha_copy, next_x_beta_geno, next_x_beta_copy;
	long cur_x_A_copy, cur_x_B_copy;
	double pi = 3.1415927; 

	*x_loglike = 0.0; 
	*place_x_loglike = 0.0;
	*y_loglike = 0.0;


	// First get the score for the X's 
	for(t = 0; t < data->seq_len-1; t++)
	{
		cur_x_alpha_geno = BIT(data->x_vec[t], 0); cur_x_alpha_copy = BITS(data->x_vec[t], 1, 3);
		cur_x_beta_geno = BIT(data->x_vec[t], 4);  cur_x_beta_copy = BITS(data->x_vec[t], 5, 7);
		next_x_alpha_geno = BIT(data->x_vec[t+1], 0); next_x_alpha_copy = BITS(data->x_vec[t+1], 1, 3);
		next_x_beta_geno = BIT(data->x_vec[t+1], 4);  next_x_beta_copy = BITS(data->x_vec[t+1], 5, 7);
		cur_x_A_copy = (1^cur_x_alpha_geno)*cur_x_alpha_copy + (1^cur_x_beta_geno)*cur_x_beta_copy;
		cur_x_B_copy = cur_x_alpha_geno*cur_x_alpha_copy + cur_x_beta_geno*cur_x_beta_copy;

		(*x_loglike) += ( log( hmm->M[cur_x_alpha_copy][next_x_alpha_copy] ) + 
				log( hmm->M[cur_x_beta_copy][next_x_beta_copy] ) );
		(place_x_loglike[t]) = 	( log( hmm->place_M[0/*cur_x_alpha_geno*/][next_x_alpha_geno][t+1] ) + 
				log( hmm->place_M[0/*cur_x_beta_geno*/][next_x_beta_geno][t+1] ) );  // The X's contribution
		(*y_loglike) -= ( (data->y_vecB[t] - hmm->MU[cur_x_A_copy][0]) * (data->y_vecB[t] - hmm->MU[cur_x_A_copy][0]) / 
			(2*hmm->SIGMA[cur_x_A_copy][0]*hmm->SIGMA[cur_x_A_copy][0]) + 
			0.5*log(2*pi) + log(hmm->SIGMA[cur_x_A_copy][0]) + 
			(data->y_vec[t] - hmm->MU[cur_x_B_copy][0]) * (data->y_vec[t] - hmm->MU[cur_x_B_copy][0]) / 
			(2*hmm->SIGMA[cur_x_B_copy][0]*hmm->SIGMA[cur_x_B_copy][0]) + 
			0.5*log(2*pi) + log(hmm->SIGMA[cur_x_B_copy][0]) );

		loglike_score += ( log( hmm->M[cur_x_alpha_copy][next_x_alpha_copy] ) + 
				log( hmm->place_M[cur_x_alpha_geno][next_x_alpha_geno][t] ) + 
				log( hmm->M[cur_x_beta_copy][next_x_beta_copy] ) + 
				log( hmm->place_M[cur_x_beta_geno][next_x_beta_geno][t] ) );  // The X's contribution


		loglike_score -= ( (data->y_vecB[t] - hmm->MU[cur_x_A_copy][0]) * (data->y_vecB[t] - hmm->MU[cur_x_A_copy][0]) / 
			(2*hmm->SIGMA[cur_x_A_copy][0]*hmm->SIGMA[cur_x_A_copy][0]) + 
			0.5*log(2*pi) + log(hmm->SIGMA[cur_x_A_copy][0]) + 
			(data->y_vec[t] - hmm->MU[cur_x_B_copy][0]) * (data->y_vec[t] - hmm->MU[cur_x_B_copy][0]) / (
			2*hmm->SIGMA[cur_x_B_copy][0]*hmm->SIGMA[cur_x_B_copy][0]) + 
			0.5*log(2*pi) + log(hmm->SIGMA[cur_x_B_copy][0]) ); // The Y's contribution
	}


	return loglike_score;

}



// Perform the Viterbi algorithm. Find the most probable path for the X's given a sequence of Y's
// Implementation is according to Rabiner's tutorial
// Everything is done with logs to avoid underflows
// Note: We added a new option which includes data from TWO output vectors !!!!! 
long Viterbi( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_X_VALS],  double *y_cond_tabsB[MAX_X_VALS],
			  long *opt_x_vec_out, double *x_loglike, double *x_place_loglike, double *y_loglike)
{

	long i, j, t;
	
	// For now take too much memory, for the special case ...
	double *lambda[MAX_X_VALS];

	// This is for debug:
//	double *lambda_x[MAX_X_VALS]; 
//	double *lambda_x_place[MAX_X_VALS];
//	double *lambda_y[MAX_X_VALS];

//	double max_x_loglike; //, max_x_place_loglike, max_y_loglike;

	long *max_indexes[MAX_X_VALS];

	double pi=3.1415927;

	double cur_log_prob, max_log_prob, cur_log_prob_lambda, max_log_prob_lambda;
	long max_index;


	*y_loglike=0; *x_loglike=0; *x_place_loglike=0; 

	// Initilize and allocate memory
	for(j = 0; j < hmm->x_dim; j++)
	{
		lambda[j] = new double[data->seq_len];
		max_indexes[j] = new long[data->seq_len];
		lambda[j][0] = log(hmm->PI[j]) + log(y_cond_tabs[j][0]);		
		max_indexes[j][0] = -1;   // maybe 0/-1 ??? Why -1 ? probably irrelevant
	}


	// Perform recursive procedure
	if(data->miss_data == 0) // no missing data
	{
		for(t = 1; t < data->seq_len; t++)
		{
			// Find the maximal lambda
			for(j = 0; j < hmm->x_dim; j++)			// state at time t
			{
				max_index = -1; max_log_prob = -9999999999999; max_log_prob_lambda = -999999999999;
				for(i = 0; i < hmm->x_dim; i++)		// state at time t-1
				{
					cur_log_prob = lambda[i][t-1] + log(hmm->M[i][j]);
					cur_log_prob_lambda = cur_log_prob +  log(y_cond_tabs[j][t]); 
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
	}
	else          // here data is missing - probably not working yet ..
		for(t = 1; t < data->seq_len; t++)
		{
			// Find the maximal lambda
			for(j = 0; j < hmm->x_dim; j++)			// state at time t
			{
				max_index = -1; max_log_prob = -9999999999999; max_log_prob_lambda = -999999999999;
				for(i = 0; i < hmm->x_dim; i++)		// state at time t-1
				{
					cur_log_prob = lambda[i][t-1] + log(hmm->M_POWERS[data->loc_diff_vec[t]][i][j]);
					cur_log_prob_lambda = cur_log_prob +  log(y_cond_tabs[j][t]); 
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


	// Termination 
	double p_max = -99999999999999.9;
	long p_max_index = -123;


	for(i = 0; i < hmm->x_dim; i++)		// state at time t-1
		if(p_max < lambda[i][t-1])
		{
			p_max = lambda[i][t-1];
			p_max_index = i;
		}


	// Find optimal path using backtracking
	opt_x_vec_out[data->seq_len-1] = p_max_index;	
	for(t = data->seq_len-2; t >= 0; t--)
		opt_x_vec_out[t] = max_indexes[opt_x_vec_out[t+1]][t+1];

	// Lastly free memory 
	for(j = 0; j < hmm->x_dim; j++)
	{
		delete lambda[j];
		delete max_indexes[j];
	}
		
	return 0; 
}






// Perform the Viterbi algorithm. Find the most probable path for the X's given a sequence of Y's
// Implementation is according to Rabiner's tutorial
// Everything is done with logs to avoid underflows
// Note: We added a new option which includes data from TWO output vectors !!!!! 
long ViterbiSNPs( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_2D_X_VALS],  double *y_cond_tabsB[MAX_2D_X_VALS],
			  long *opt_x_vec_out) 
{
	long i, j, t, i1,i2,j1,j2, i_geno, j_geno, total_i_index, total_j_index;
	long i_geno1, i_geno2, j_geno1, j_geno2, total_i_index0, total_i_index1;
	
	long t_xy, t_xz, t_yz; // for argmax finding

	// Convention: First I, then J. First Copy, then Genotype. First 1, then 2.
	// These temporary tables are for faster inplementation (instead of summing over all I's and J's
	double K_tab1[3][3][2] [2];
	double K_tab2[3][3] [2][2];
	double K_tab3[3] [3][2][2];

	long Ind_tab1[3][3][2] [2];
	long Ind_tab2[3][3] [2][2];
	long Ind_tab3[3] [3][2][2];
	
	// For now take too much memory, for the special case ...
	double *lambda[36];


	long *max_indexes[36];

	double pi=3.1415927;

	double cur_log_prob, max_log_prob, cur_log_prob_lambda, max_log_prob_lambda;
	long max_index;



	// Initilize and allocate memory
	for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	 // two genotypes indexes
		for(j1 = 0; j1 < hmm->x_dim; j1++)			// two copy numbers
			for(j2 = 0; j2 < hmm->x_dim; j2++)
			{
				total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
				lambda[total_j_index] = new double[data->seq_len];
				max_indexes[total_j_index] = new long[data->seq_len];
				{
					if(hmm->gauss_dim == 1) // old version
						lambda[total_j_index][0] = log(hmm->PI[j1]) + log(hmm->PI[j2]) + 
							log(y_cond_tabs[A_copy_tab[total_j_index]][0]) + 
							log(y_cond_tabsB[B_copy_tab[total_j_index]][0]);								
					else // new 2d version
						lambda[total_j_index][0] = log(hmm->PI[j1]) + log(hmm->PI[j2]) + 
							log(y_cond_tabs[A_copy_tab[total_j_index]+5*B_copy_tab[total_j_index]][0]);				
					max_indexes[total_j_index][0] = -1;   // maybe 0/-1 ???
				}
			}

	// Perform recursive procedure
	if(data->miss_data == 0) // no missing data
	{
		for(t = 1; t < data->seq_len; t++)
		{
/******/
			// Find the maximal lambda
//////////////////////////////////// PHASE 1 ////////////////////////////
			for(j_geno2 = 0; j_geno2 < hmm->x_dim2; j_geno2++)	 // genotype 2 indexes
				for(i1 = 0; i1 < hmm->x_dim; i1++)	
					for(i2 = 0; i2 < hmm->x_dim; i2++)	
						for(i_geno1 = 0; i_geno1 < hmm->x_dim2; i_geno1++)	
						{
							total_i_index0 = multi_dim_total_index_tab[i1][i_geno1][i2][0]; 
							total_i_index1 = multi_dim_total_index_tab[i1][i_geno1][i2][1]; 
							K_tab1[i1][i2][i_geno1][j_geno2] = MAX( lambda[total_i_index0][t-1] + 
								log(hmm->place_M[0][j_geno2][t-1]),  
							lambda[total_i_index1][t-1] + 
								log(hmm->place_M[1][j_geno2][t-1]) ); // max over i_geno2

							Ind_tab1[i1][i2][i_geno1][j_geno2] = ARGMAX( lambda[total_i_index0][t-1] + 
								log(hmm->place_M[0][j_geno2][t-1]),  
							lambda[total_i_index1][t-1] + 
								log(hmm->place_M[1][j_geno2][t-1]) ); // argmax over i_geno2

						}
//////////////////////////////////// PHASE 2 ////////////////////////////
			for(j_geno1 = 0; j_geno1 < hmm->x_dim2; j_geno1++)	 // genotype 1 indexes
				for(j_geno2 = 0; j_geno2 < hmm->x_dim2; j_geno2++)	 // genotype 2 indexes
					for(i1 = 0; i1 < hmm->x_dim; i1++)	
						for(i2 = 0; i2 < hmm->x_dim; i2++)	
						{
							K_tab2[i1][i2][j_geno1][j_geno2] = MAX( K_tab1[i1][i2][0][j_geno2] + 
								log(hmm->place_M[0][j_geno1][t-1]),  
							K_tab1[i1][i2][1][j_geno2]  + 
								log(hmm->place_M[1][j_geno1][t-1]) ); // max over i_geno1

							Ind_tab2[i1][i2][j_geno1][j_geno2] = ARGMAX( K_tab1[i1][i2][0][j_geno2] + 
								log(hmm->place_M[0][j_geno1][t-1]),  
							K_tab1[i1][i2][1][j_geno2]  + 
								log(hmm->place_M[1][j_geno1][t-1]) ); // argmax over i_geno1
						}
//////////////////////////////////// PHASE 3 ////////////////////////////
			for(j_geno1 = 0; j_geno1 < hmm->x_dim2; j_geno1++)	 // genotype 1 indexes
				for(j_geno2 = 0; j_geno2 < hmm->x_dim2; j_geno2++)	 // genotype 2 indexes
					for(i1 = 0; i1 < hmm->x_dim; i1++)	
						for(j2 = 0; j2 < hmm->x_dim; j2++)	
						{
							K_tab3[i1][j2][j_geno1][j_geno2] = MAX( MAX( K_tab2[i1][0][j_geno1][j_geno2] + 
								log(hmm->M[0][j2]),  
							K_tab2[i1][1][j_geno1][j_geno2] + 
								log(hmm->M[1][j2])), 
							K_tab2[i1][2][j_geno1][j_geno2] + 
								log(hmm->M[2][j2]) ); // max over i2

							ARGMAX3( K_tab2[i1][0][j_geno1][j_geno2] + log(hmm->M[0][j2]),  
									 K_tab2[i1][1][j_geno1][j_geno2] + log(hmm->M[1][j2]), 
									 K_tab2[i1][2][j_geno1][j_geno2] + log(hmm->M[2][j2]), 
									 Ind_tab3[i1][j2][j_geno1][j_geno2] ); // argmax over i2

///							Ind_tab3[i1][j2][j_geno1][j_geno2] = ARGMAX(110,10);
							
///							ARGMAX3(2,1,0, Ind_tab3[i1][j2][j_geno1][j_geno2] );
						}
//////////////////////////////////// PHASE 4 ////////////////////////////
			for(j_geno1 = 0; j_geno1 < hmm->x_dim2; j_geno1++)	 // genotype 1 indexes
				for(j_geno2 = 0; j_geno2 < hmm->x_dim2; j_geno2++)	 // genotype 2 indexes
					for(j1 = 0; j1 < hmm->x_dim; j1++)	
						for(j2 = 0; j2 < hmm->x_dim; j2++)	
						{
							total_j_index = multi_dim_total_index_tab[j1][j_geno1][j2][j_geno2]; 
							if(hmm->gauss_dim == 1) // old version
								lambda[total_j_index][t] = ( MAX( MAX( K_tab3[0][j2][j_geno1][j_geno2] + 
									log(hmm->M[0][j1]), 
								K_tab3[1][j2][j_geno1][j_geno2] +
									log(hmm->M[1][j1])), 
								K_tab3[2][j2][j_geno1][j_geno2] + 
									log(hmm->M[2][j1])) ) +  // max over i1
										log(y_cond_tabs[A_copy_tab[total_j_index]][t]) + 
										log(y_cond_tabsB[B_copy_tab[total_j_index]][t]);	
							else
								lambda[total_j_index][t] = ( MAX( MAX( K_tab3[0][j2][j_geno1][j_geno2] + 
									log(hmm->M[0][j1]), 
								K_tab3[1][j2][j_geno1][j_geno2] +
									log(hmm->M[1][j1])), 
								K_tab3[2][j2][j_geno1][j_geno2] + 
									log(hmm->M[2][j1])) ) +  // max over i1
										log(y_cond_tabs[A_copy_tab[total_j_index]+5*B_copy_tab[total_j_index]][t]); 
							
							ARGMAX3( K_tab3[0][j2][j_geno1][j_geno2] + log(hmm->M[0][j1]), 
									 K_tab3[1][j2][j_geno1][j_geno2] + log(hmm->M[1][j1]), 
									 K_tab3[2][j2][j_geno1][j_geno2] + log(hmm->M[2][j1]), 
									 i1 ); // argmax over i1
							i2 = Ind_tab3[i1][j2][j_geno1][j_geno2];
							i_geno1 = Ind_tab2[i1][i2][j_geno1][j_geno2];
							i_geno2 = Ind_tab1[i1][i2][i_geno1][j_geno2];
							max_indexes[total_j_index][t] = multi_dim_total_index_tab[i1][i_geno1][i2][i_geno2];   // take the maximal index. How the hell do we know it?
						}
////////////////////////////////////////////////////////////////
		}			// end loop on t
	}
	else          // here data is missing - Certainly not working yet, no ENERGY !!!
		for(t = 1; t < data->seq_len; t++)
		{
			// Find the maximal lambda
			for(j = 0; j < hmm->x_dim; j++)			// state at time t
			{
				max_index = -1; max_log_prob = -9999999999999; max_log_prob_lambda = -999999999999;
				for(i = 0; i < hmm->x_dim; i++)		// state at time t-1
				{
					cur_log_prob = lambda[i][t-1] + log(hmm->M_POWERS[data->loc_diff_vec[t]][i][j]);
					cur_log_prob_lambda = cur_log_prob +  log(y_cond_tabs[j][t]); 
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


	// Termination 
	double p_max = -99999999999999.9;
	long p_max_index = -123;


	for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)		// state at time t-1
		for(i1 = 0; i1 < hmm->x_dim; i1++)		// state at time t-1
			for(i2 = 0; i2 < hmm->x_dim; i2++)		// state at time t-1
			{
				total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 
				if(p_max < lambda[total_i_index][t-1])
				{
					p_max = lambda[total_i_index][t-1];
					p_max_index = total_i_index;
				}
			}

	// Find optimal path using backtracking
	opt_x_vec_out[data->seq_len-1] = p_max_index;	
	for(t = data->seq_len-2; t >= 0; t--)
		opt_x_vec_out[t] = max_indexes[opt_x_vec_out[t+1]][t+1];

///	return 111;	

	// Lastly free memory 
	for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	 // two genotypes indexes
		for(j1 = 0; j1 < hmm->x_dim; j1++)			// two copy numbers
			for(j2 = 0; j2 < hmm->x_dim; j2++)				
			{
				total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)];
				delete lambda[total_j_index];
				delete max_indexes[total_j_index];
			}
		
	return 0; 
}



// Perform the Forward algorithm. Find the alpha coefficients and the prob. of Y (summing over all possible X's)
// Implementation is according to Rabiner's tutorial.
// We need to do 'scaling' of the alphas to avoid underflows ..
// alpha[j][t] = Pr(x_t = j | y_1,..,y_t)
long forward( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_X_VALS], 
			  double *alpha[MAX_X_VALS], double *scale, double *scale_cum,  double *scale_exp, double *y_vec_log_prob_out)
{
	long i, j, t;
	

	// Initilize
	for(j = 0; j < hmm->x_dim; j++)
		alpha[j][0] = hmm->PI[j] * y_cond_tabs[j][0];    // use the already computed b tables 
	scale[0] = 0; scale_exp[0] = 1;


	// Perform recursive procedure	
	if(data->miss_data == 0)  // No missing data - the simple case 
		for(t = 1; t < data->seq_len; t++)
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
					alpha[j][t] *= y_cond_tabs[j][t];
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
					alpha[j][t] *= y_cond_tabs[j][t];
				}
			}
		}
	else   // here we must consider missing data !!! 
		for(t = 1; t < data->seq_len; t++)
		{
			// Compute the next alphas. Note : We need scaling here !!!
			if((t%DO_SCALE)==DO_SCALE-1)
			{
				scale_exp[t]=0;
				for(j = 0; j < hmm->x_dim; j++)			// state at time t
				{
					alpha[j][t] = alpha[0][t-1] * hmm->M_POWERS[data->loc_diff_vec[t]][0][j]; // take the 1st one ..	
					for(i = 1; i < hmm->x_dim; i++)		// state at time t-1
						alpha[j][t] += alpha[i][t-1] * hmm->M_POWERS[data->loc_diff_vec[t]][i][j];
					alpha[j][t] *= y_cond_tabs[j][t];
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
						alpha[j][t] = alpha[0][t-1] * hmm->M_POWERS[data->loc_diff_vec[t]][0][j]; // take the 1st one ..	
		
					for(i = 1; i < hmm->x_dim; i++)		// state at time t-1
						alpha[j][t] += alpha[i][t-1] * hmm->M_POWERS[data->loc_diff_vec[t]][i][j];
					alpha[j][t] *= y_cond_tabs[j][t];
				}
			}
		}


	// Get the scaling of all alpha's : 
	scale_cum[0] = scale[0];
	for(t = 1; t < data->seq_len; t++)
		scale_cum[t] = scale_cum[t-1] + scale[t];

	*y_vec_log_prob_out = scale_cum[data->seq_len-1]; // output the total scaling in log

	return 0; 
}






// New: Try to save time AND get the same results ...
// A special function for performing the Forward algorithm in the complicated SNPs model. 
// Find the alpha coefficients and the prob. of Y (summing over all possible X's)
// Implementation is according to Rabiner's tutorial.
// We need to do 'scaling' of the alphas to avoid underflows ..
// The inputs are: 
// hmm - the model
// data - the data
// y_cond_tabs - pre-calculation of y's conditional probabilitiyes
// y_cond_tabsB - pre-calculation of y's B conditional probabilitiyes
//
// The outputs are:
// alpha - the local conditional probabilites
// scale - scaling parameters for the alphas.
// scale_cum - cumulative sum of scale
// scale_exp - exponent of scale
// y_vec_log_prob_out - log probabilities of the y's 
long forwardSNPs( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS], 
			  double *alpha[36], double *scale, double *scale_cum,  double *scale_exp, double *y_vec_log_prob_out)
{
	long i, j, t, i1, j1, i2, j2, i_geno, j_geno, total_i_index, total_j_index;
	long i_geno1, j_geno1, j_geno2, total_i_index0, total_i_index1;
	
	// Convention: First I, then J. First Copy, then Genotype. First 1, then 2.
	// These temporary tables are for faster inplementation (instead of summing over all I's and J's
	double K_tab1[3][3][2] [2];
	double K_tab2[3][3] [2][2];
	double K_tab3[3] [3][2][2];


	// Initilize
	for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	 // two genotypes indexes
		for(j1 = 0; j1 < hmm->x_dim; j1++)			// two copy numbers
			for(j2 = 0; j2 < hmm->x_dim; j2++)
			{
				total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
				// use the already computed b tables. Compensate for the possible 4 different genotypes at the beginning
				if(hmm->gauss_dim == 1) // old version
					alpha[total_j_index][0] = hmm->PI[j1] * y_cond_tabs[A_copy_tab[total_j_index]][0] * 
											  hmm->PI[j2] * y_cond_tabsB[B_copy_tab[total_j_index]][0] * 0.25;  
				else // new 2d version
					alpha[total_j_index][0] = hmm->PI[j1] * y_cond_tabs[A_copy_tab[total_j_index]+5*B_copy_tab[total_j_index]][0]; 

			}
	scale[0] = 0; scale_exp[0] = 1;


	// Perform recursive procedure	
	if(data->miss_data == 0)  // No missing data - the simple case 
	{
		for(t = 1; t < data->seq_len; t++)
		{
			// Compute the next alphas. Note : We need scaling here !!!
			if((t%DO_SCALE)==DO_SCALE-1)
			{
				scale_exp[t]=0;
//////////////////////////////////// PHASE 1 ////////////////////////////
				for(j_geno2 = 0; j_geno2 < hmm->x_dim2; j_geno2++)	 // genotype 2 indexes
					for(i1 = 0; i1 < hmm->x_dim; i1++)	
						for(i2 = 0; i2 < hmm->x_dim; i2++)	
							for(i_geno1 = 0; i_geno1 < hmm->x_dim2; i_geno1++)	
							{
								total_i_index0 = multi_dim_total_index_tab[i1][i_geno1][i2][0]; 
								total_i_index1 = multi_dim_total_index_tab[i1][i_geno1][i2][1]; 
								K_tab1[i1][i2][i_geno1][j_geno2] = alpha[total_i_index0][t-1] * 
									hmm->place_M[0][j_geno2][t-1] + 
								alpha[total_i_index1][t-1] * 
									hmm->place_M[1][j_geno2][t-1]; // sum over i_geno2
							}

//////////////////////////////////// PHASE 2 ////////////////////////////
				for(j_geno1 = 0; j_geno1 < hmm->x_dim2; j_geno1++)	 // genotype 1 indexes
					for(j_geno2 = 0; j_geno2 < hmm->x_dim2; j_geno2++)	 // genotype 2 indexes
						for(i1 = 0; i1 < hmm->x_dim; i1++)	
							for(i2 = 0; i2 < hmm->x_dim; i2++)	
							{
								K_tab2[i1][i2][j_geno1][j_geno2] = K_tab1[i1][i2][0][j_geno2] * 
									hmm->place_M[0][j_geno1][t-1] + 
								K_tab1[i1][i2][1][j_geno2]  * 
									hmm->place_M[1][j_geno1][t-1]; // sum over i_geno1
							}
//////////////////////////////////// PHASE 3 ////////////////////////////
				for(j_geno1 = 0; j_geno1 < hmm->x_dim2; j_geno1++)	 // genotype 1 indexes
					for(j_geno2 = 0; j_geno2 < hmm->x_dim2; j_geno2++)	 // genotype 2 indexes
						for(i1 = 0; i1 < hmm->x_dim; i1++)	
							for(j2 = 0; j2 < hmm->x_dim; j2++)	
							{
								K_tab3[i1][j2][j_geno1][j_geno2] = K_tab2[i1][0][j_geno1][j_geno2] * 
									hmm->M[0][j2] + 
								K_tab2[i1][1][j_geno1][j_geno2] * 
									hmm->M[1][j2] +
								K_tab2[i1][2][j_geno1][j_geno2] * 
									hmm->M[2][j2]; // sum over i2
							}
//////////////////////////////////// PHASE 4 ////////////////////////////
				for(j_geno1 = 0; j_geno1 < hmm->x_dim2; j_geno1++)	 // genotype 1 indexes
					for(j_geno2 = 0; j_geno2 < hmm->x_dim2; j_geno2++)	 // genotype 2 indexes
						for(j1 = 0; j1 < hmm->x_dim; j1++)	
							for(j2 = 0; j2 < hmm->x_dim; j2++)	
							{
								total_j_index = multi_dim_total_index_tab[j1][j_geno1][j2][j_geno2]; 
								alpha[total_j_index][t] = (K_tab3[0][j2][j_geno1][j_geno2] * 
									hmm->M[0][j1] + 
								K_tab3[1][j2][j_geno1][j_geno2] * 
									hmm->M[1][j1] +
								K_tab3[2][j2][j_geno1][j_geno2] * 
									hmm->M[2][j1]);
								
								if(hmm->gauss_dim == 1) // old version
									alpha[total_j_index][t] *=  
										( y_cond_tabs[A_copy_tab[total_j_index]][t] * 
										y_cond_tabsB[B_copy_tab[total_j_index]][t] );
								else // here gaussian 2-dim
									alpha[total_j_index][t] *=  
										y_cond_tabs[A_copy_tab[total_j_index]+5*B_copy_tab[total_j_index]][t];
								scale_exp[t] += alpha[total_j_index][t]; 	
							}
////////////////////////////////////////////////////////////////
						
				// Now do scaling	
				for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	 // two genotypes indexes
					for(j1 = 0; j1 < hmm->x_dim; j1++)			// two copy numbers
						for(j2 = 0; j2 < hmm->x_dim; j2++)					
						{
							total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
							alpha[total_j_index][t] /= scale_exp[t];
						}
				scale[t] = log(scale_exp[t]); // transfer to logs ..
			}
			else  // No Scaling !! does it work? unknown .. 
			{
				scale[t] = 0; scale_exp[t] = 1;
				for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	 // two genotypes indexes
					for(j1 = 0; j1 < hmm->x_dim; j1++)			// two copy numbers
						for(j2 = 0; j2 < hmm->x_dim; j2++)					
						{
							total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
							alpha[total_j_index][t] = 0; // take the 1st one ..	
							for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)	
								for(i1 = 0; i1 < hmm->x_dim; i1++)	
									for(i2 = 0; i2 < hmm->x_dim; i2++)	
									{
										total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 
										alpha[total_j_index][t] += alpha[total_i_index][t-1] * 
											hmm->M[i1][j1] * hmm->M[i2][j2] * 
											hmm->place_M[BIT(i_geno,0)][BIT(j_geno,0)][t-1] * 
											hmm->place_M[BIT(i_geno,1)][BIT(j_geno,1)][t-1];
									}
							if(hmm->gauss_dim == 1) // old version
								alpha[total_j_index][t] *= ( y_cond_tabs[A_copy_tab[total_j_index]][t] * 
															 y_cond_tabsB[B_copy_tab[total_j_index]][t] );			
							else
								alpha[total_j_index][t] *= y_cond_tabs[A_copy_tab[total_j_index]+5*B_copy_tab[total_j_index]][t];
						}
			}
		}
	}
	else   // here we must consider missing data. Not working for now due to lack of ENERGY !!! !!! 
		for(t = 1; t < data->seq_len; t++)
		{
			// Compute the next alphas. Note : We need scaling here !!!
			if((t%DO_SCALE)==DO_SCALE-1)
			{
				scale_exp[t]=0;
				for(j = 0; j < hmm->x_dim; j++)			// state at time t
				{
					alpha[j][t] = alpha[0][t-1] * hmm->M_POWERS[data->loc_diff_vec[t]][0][j]; // take the 1st one ..	
					for(i = 1; i < hmm->x_dim; i++)		// state at time t-1
						alpha[j][t] += alpha[i][t-1] * hmm->M_POWERS[data->loc_diff_vec[t]][i][j];
					alpha[j][t] *= y_cond_tabs[j][t];
					scale_exp[t] += alpha[j][t];
				}
		
				// Now do scaling	
				for(j = 0; j < hmm->x_dim; j++)
					alpha[j][t] /= scale_exp[t];
				scale[t] = log(scale_exp[t]); // transfer to logs ..
			}
			else  // No Scaling !!
			{
				scale[t] = 0; scale_exp[t] = 1;
				for(j = 0; j < hmm->x_dim; j++)			// state at time t
				{
					alpha[j][t] = alpha[0][t-1] * hmm->M_POWERS[data->loc_diff_vec[t]][0][j]; // take the 1st one ..	
					for(i = 1; i < hmm->x_dim; i++)		// state at time t-1
						alpha[j][t] += alpha[i][t-1] * hmm->M_POWERS[data->loc_diff_vec[t]][i][j];
					alpha[j][t] *= y_cond_tabs[j][t];
				}
			}
		}  // finished if(missing_data) part


	// Get the scaling of all alpha's : 
	scale_cum[0] = scale[0];
	for(t = 1; t < data->seq_len; t++)
		scale_cum[t] = scale_cum[t-1] + scale[t];

	*y_vec_log_prob_out = scale_cum[data->seq_len-1]; // output the total scaling in log

	return 0; 
}





// Old: get the same results in a 'standard' way - just good for 'pedagogical reasons'
// A special function for performing the Forward algorithm in the complicated SNPs model. 
// Find the alpha coefficients and the prob. of Y (summing over all possible X's)
// Implementation is according to Rabiner's tutorial.
// We need to do 'scaling' of the alphas to avoid underflows ..
long SLOWforwardSNPs( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS], 
			  double *alpha[36], double *scale, double *scale_cum,  double *scale_exp, double *y_vec_log_prob_out)
{
	long i, j, t, i1, j1, i2, j2, i_geno, j_geno, total_i_index, total_j_index;
	long i_geno1, j_geno1, j_geno2, total_i_index0, total_i_index1;
	


	// Initilize
	for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	 // two genotypes indexes
		for(j1 = 0; j1 < hmm->x_dim; j1++)			// two copy numbers
			for(j2 = 0; j2 < hmm->x_dim; j2++)
			{
				total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
				// use the already computed b tables. Compensate for the possible 4 different genotypes at the beginning
				alpha[total_j_index][0] = hmm->PI[j1] * y_cond_tabs[A_copy_tab[total_j_index]][0] * 
										  hmm->PI[j2] * y_cond_tabsB[B_copy_tab[total_j_index]][0] * 0.25;  
			}
	scale[0] = 0; scale_exp[0] = 1;


	// Perform recursive procedure	
	if(data->miss_data == 0)  // No missing data - the simple case 
	{
		for(t = 1; t < data->seq_len; t++)
		{
			// Compute the next alphas. Note : We need scaling here !!!
			if((t%DO_SCALE)==DO_SCALE-1)
			{
				scale_exp[t]=0;
				for(j_geno = 0; j_geno < hmm->x_dim2*hmm->x_dim2; j_geno++)	 // genotype 2 indexes
					for(j1 = 0; j1 < hmm->x_dim; j1++)	
						for(j2 = 0; j2 < hmm->x_dim; j2++)	
						{
							total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
							alpha[total_j_index][t] = 0;
							for(i_geno = 0; i_geno < hmm->x_dim2*hmm->x_dim2; i_geno++)	 // genotype 2 indexes
								for(i1 = 0; i1 < hmm->x_dim; i1++)	
									for(i2 = 0; i2 < hmm->x_dim; i2++)	
									{
										total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 
										alpha[total_j_index][t] += ( alpha[total_i_index][t-1] * 
											hmm->M[i1][j1] * hmm->M[i2][j2] * 
											hmm->place_M[BIT(i_geno,0)][BIT(j_geno,0)][t-1] * 
											hmm->place_M[BIT(i_geno,1)][BIT(j_geno,1)][t-1] );
									}

							alpha[total_j_index][t] *=
								( y_cond_tabs[A_copy_tab[total_j_index]][t] * 
								y_cond_tabsB[B_copy_tab[total_j_index]][t] );	
							scale_exp[t] += alpha[total_j_index][t]; 	
						}
						
				// Now do scaling	
				for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	 // two genotypes indexes
					for(j1 = 0; j1 < hmm->x_dim; j1++)			// two copy numbers
						for(j2 = 0; j2 < hmm->x_dim; j2++)					
						{
							total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
							alpha[total_j_index][t] /= scale_exp[t];
						}
				scale[t] = log(scale_exp[t]); // transfer to logs ..
			}
		}
	}

	// Get the scaling of all alpha's : 
	scale_cum[0] = scale[0];
	for(t = 1; t < data->seq_len; t++)
		scale_cum[t] = scale_cum[t-1] + scale[t];

	*y_vec_log_prob_out = scale_cum[data->seq_len-1]; // output the total scaling in log

	return 0; 
}









// Perform the Backward algorithm. Find the beta coefficients and the prob. of Y (summing over all possible X's)
// Implementation is according to Rabiner's tutorial
// Note : We assume here that the 'scaling' coefficients are already determined by the forward algorithm
long backward( hmm_model *hmm, hmm_data *data, double *scale_exp, double *y_cond_tabs[MAX_X_VALS],
			  double *beta[MAX_X_VALS])
{
	long i, j, t;
	

	// Initilize
	for(i = 0; i < hmm->x_dim; i++)
		beta[i][data->seq_len-1] = 1 / scale_exp[data->seq_len-1]; // Here we start with one scaled.		
	
	if(data->miss_data == 0)  // No missing data
		for(t = data->seq_len-2; t >= 0; t--)      // Compute the next betas. Note : We need scaling here !!!
			for(i = 0; i < hmm->x_dim; i++)			// state at time t
			{
				beta[i][t] = beta[0][t+1] * hmm->M[i][0] *  y_cond_tabs[0][t+1]; // take the 1st one ..
				for(j = 1; j < hmm->x_dim; j++)		// state at time t+1
					beta[i][t] += beta[j][t+1] * hmm->M[i][j] * y_cond_tabs[j][t+1]; 
				beta[i][t] /= scale_exp[t];  // already apply scaling !!!
			}
	else // Consider missing data
		for(t = data->seq_len-2; t >= 0; t--)      // Compute the next betas. Note : We need scaling here !!!
			for(i = 0; i < hmm->x_dim; i++)			// state at time t
			{
				beta[i][t] = beta[0][t+1] * hmm->M_POWERS[data->loc_diff_vec[t+1]][i][0] *  
					y_cond_tabs[0][t+1]; // take the 1st one ..
				for(j = 1; j < hmm->x_dim; j++)		// state at time t+1
					beta[i][t] += beta[j][t+1] * hmm->M_POWERS[data->loc_diff_vec[t+1]][i][j] * 
					y_cond_tabs[j][t+1]; 
				beta[i][t] /= scale_exp[t];  // already apply scaling !!!
			}

	// Note : The M matrices used here must be the same as the matrices used in the forward algorithm,
	// otherwise we have a serious scaling problem ..

	return 0; 

}




// A special function for performing the Backward algorithm in the complicated SNPs model. 
// Find the beta coefficients and the prob. of Y (summing over all possible X's)
// Implementation is according to Rabiner's tutorial
// Note : We assume here that the 'scaling' coefficients are already determined by the forward algorithm
// The output is the beta vec/matrix
long backwardSNPs( hmm_model *hmm, hmm_data *data, double *scale_exp, 
				   double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS],
				   double *beta[36])
{
	long i, j, t, i1, i2, j1, j2, total_i_index;
	long i_geno1, i_geno2, j_geno1, total_j_index0, total_j_index1;

	// Convention: First Copy, then Genotype. First 1, then 2. First I, then J. 
	// These temporary tables are for faster inplementation (instead of summing over all I's and J's
	double K_tab1[3][3][2] [2];
	double K_tab2[3][3] [2][2];
	double K_tab3[3] [3][2][2];

	

	// Initilize
	for(i=0; i < (hmm->x_dim*hmm->x_dim*hmm->x_dim2*hmm->x_dim2); i++)
		beta[i][data->seq_len-1] = 1 / scale_exp[data->seq_len-1]; // Here we start with one scaled.

	if(data->miss_data == 0)  // No missing data
		for(t = data->seq_len-2; t >= 0; t--)      // Compute the next betas. Note : We need scaling here !!!
		{
//////////////////////////////////// PHASE 1 ////////////////////////////
			for(j_geno1 = 0; j_geno1 < hmm->x_dim2; j_geno1++)	 // genotype 2 indexes
				for(j1 = 0; j1 < hmm->x_dim; j1++)	
					for(j2 = 0; j2 < hmm->x_dim; j2++)	
						for(i_geno2 = 0; i_geno2 < hmm->x_dim2; i_geno2++)	
						{
							total_j_index0 = multi_dim_total_index_tab[j1][j_geno1][j2][0]; 
							total_j_index1 = multi_dim_total_index_tab[j1][j_geno1][j2][1]; 

							if(hmm->gauss_dim == 1) // standard old way ... 
								K_tab1[j1][j2][j_geno1][i_geno2] = beta[total_j_index0][t+1] * 
									y_cond_tabs[A_copy_tab[total_j_index0]][t+1] * 
									y_cond_tabsB[B_copy_tab[total_j_index0]][t+1] *
										hmm->place_M[i_geno2][0][t] + 
										beta[total_j_index1][t+1] * 
									y_cond_tabs[A_copy_tab[total_j_index1]][t+1] * 
									y_cond_tabsB[B_copy_tab[total_j_index1]][t+1] *
										hmm->place_M[i_geno2][1][t]; // sum over j_geno2
							else // new 2-d gaussians
								K_tab1[j1][j2][j_geno1][i_geno2] = beta[total_j_index0][t+1] * 
									y_cond_tabs[A_copy_tab[total_j_index0]+5*B_copy_tab[total_j_index0]][t+1] * 
										hmm->place_M[i_geno2][0][t] + 
										beta[total_j_index1][t+1] * 
									y_cond_tabs[A_copy_tab[total_j_index1]+5*B_copy_tab[total_j_index1]][t+1] * 
										hmm->place_M[i_geno2][1][t]; // sum over j_geno2

						}
//////////////////////////////////// PHASE 2 ////////////////////////////
			for(i_geno1 = 0; i_geno1 < hmm->x_dim2; i_geno1++)	 // genotype 1 indexes
				for(j1 = 0; j1 < hmm->x_dim; j1++)	
					for(j2 = 0; j2 < hmm->x_dim; j2++)	
						for(i_geno2 = 0; i_geno2 < hmm->x_dim2; i_geno2++)	
						{
							K_tab2[j1][j2][i_geno1][i_geno2] = K_tab1[j1][j2][0][i_geno2] * 
								hmm->place_M[i_geno1][0][t] + 
							K_tab1[j1][j2][1][i_geno2]  * 
								hmm->place_M[i_geno1][1][t]; // sum over j_geno1
						}
//////////////////////////////////// PHASE 3 ////////////////////////////
			for(i_geno1 = 0; i_geno1 < hmm->x_dim2; i_geno1++)	 // genotype 1 indexes
				for(j1 = 0; j1 < hmm->x_dim; j1++)	
					for(i2 = 0; i2 < hmm->x_dim; i2++)	
						for(i_geno2 = 0; i_geno2 < hmm->x_dim2; i_geno2++)	
						{
							K_tab3[j1][i2][i_geno1][i_geno2] = K_tab2[j1][0][i_geno1][i_geno2] * 
								hmm->M[i2][0] + 
							K_tab2[j1][1][i_geno1][i_geno2] * 
								hmm->M[i2][1] +
							K_tab2[j1][2][i_geno1][i_geno2] * 
								hmm->M[i2][2]; // sum over j2
						}
//////////////////////////////////// PHASE 4 ////////////////////////////
			for(i_geno1 = 0; i_geno1 < hmm->x_dim2; i_geno1++)	 // genotype 1 indexes
				for(i1 = 0; i1 < hmm->x_dim; i1++)	
					for(i2 = 0; i2 < hmm->x_dim; i2++)	
						for(i_geno2 = 0; i_geno2 < hmm->x_dim2; i_geno2++)	
						{
							total_i_index = multi_dim_total_index_tab[i1][i_geno1][i2][i_geno2]; 
							beta[total_i_index][t] = (K_tab3[0][i2][i_geno1][i_geno2] * 
								hmm->M[i1][0] + 
							K_tab3[1][i2][i_geno1][i_geno2] * 
								hmm->M[i1][1] +
							K_tab3[2][i2][i_geno1][i_geno2] * 
								hmm->M[i1][2]);  // sum over j1
							beta[total_i_index][t] /= scale_exp[t];  // already apply scaling !!!	
						}

		}

	else // Consider missing data. The part below is bad and not working. We didn't have ENERGY to change it !!!
		for(t = data->seq_len-2; t >= 0; t--)      // Compute the next betas. Note : We need scaling here !!!
			for(i = 0; i < hmm->x_dim; i++)			// state at time t
			{
				beta[i][t] = beta[0][t+1] * hmm->M_POWERS[data->loc_diff_vec[t+1]][i][0] *  
					y_cond_tabs[0][t+1]; // take the 1st one ..
				for(j = 1; j < hmm->x_dim; j++)		// state at time t+1
					beta[i][t] += beta[j][t+1] * hmm->M_POWERS[data->loc_diff_vec[t+1]][i][j] * 
					y_cond_tabs[j][t+1]; 
				beta[i][t] /= scale_exp[t];  // already apply scaling !!!
			}

	// Note : The M matrices used here must be the same as the matrices used in the forward algorithm,
	// otherwise we have a serious scaling problem ..

	return 0; 
}



// Here we find the marginal prob. of each X variable, This is an alternative to the Viterbi algorithm,
// Note that if interpreted incorrectly it can give 'illegal path' !
// We assume here we already done the forward-backward procedure !!
// The meaning of the output variables is the following:
// gamma[i][t] = Pr(X_t == i | Y_1,..,Y_n)
// phi[i][j][t] = Pr(X_t == i, Mixture_Y_t = j | Y_1,...,Y_n)
long ComputeMarginalGamma(hmm_model *hmm,  hmm_data *data, double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS],
						 double *alpha[MAX_X_VALS], double *beta[MAX_X_VALS], double *scale_exp, double *y_cond_tabs[MAX_X_VALS],
						 double *gamma[MAX_X_VALS], double *phi[MAX_X_VALS][MAX_Y_VALS])
{
	long i, j, t;

	for(i = 0; i < hmm->x_dim; i++)
		for(t = 0; t < data->seq_len; t++)
			gamma[i][t] = alpha[i][t] * beta[i][t] * scale_exp[t]; // here scale is already in exp !!! // y_vec_log_prob

	// here we compute gamma gamma !!  we use normal density here !!
	if(hmm->y_type == CONTINUOUS)
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				for(t = 0; t < data->seq_len; t++)
					phi[i][j][t] = gamma[i][t] * ///// ALREADY hmm->N[i][j] *  ONE_OVER_2SQRT_PI * 
						y_mu_square_exp_vecs[i][j][t] / y_cond_tabs[i][t]; 	
	return 0;
}


// A special function for finding the marginal prob. of each X variable, in the complicated SNPs HMM.
// This is an alternative to the Viterbi algorithm.
// Note that if interpreted incorrectly it can give 'illegal path' !
// We assume here we already done the 'special' forward-backward procedure !!
// Variables: 
// gamma[i][t] = Pr(X_t=i | Y_1,..,Y_N)
// phi[i][j][t] = Pr(X_t=i, Mixture_Y_t=j | Y_1,..,Y_N) (irrelevant for DIM_Y = 1)
long ComputeMarginalGammaSNPs(hmm_model *hmm,  hmm_data *data, 
						 double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS], double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS],
						 double *alpha[36], double *beta[36], double *scale_exp, 
						 double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS],
						 double *gamma[36], double *phi[36][MAX_Y_VALS])
{
	long t;
	long total_i_index;
	// compute only the relevant meaningful possibilities, and not the whole 36, as most of them are meaningless.

	
	for(total_i_index = 0; total_i_index < 36; total_i_index++)		// state at time t-1
		for(t = 0; t < data->seq_len; t++)
			gamma[total_i_index][t] = alpha[total_i_index][t] * beta[total_i_index][t] * scale_exp[t]; // here scale is already in exp !!! 


	// here we compute gamma gamma (phi) !!  we do not know yet what to do here !!
	if(hmm->y_type == CONTINUOUS)
		for(total_i_index = 0; total_i_index < 36; total_i_index++)		// state at time t-1
		{
#ifdef Y_DIM
			if(hmm->gauss_dim == 1) // standard old 
				for(j = 0; j < hmm->y_dim; j++)
					for(t = 0; t < data->seq_len; t++)
						phi[total_i_index][j][t] = gamma[total_i_index][t] *  
						(y_mu_square_exp_vecs[A_copy_tab[total_i_index]][j][t] / y_cond_tabs[A_copy_tab[total_i_index]][t]) * 
						(y_mu_square_exp_vecsB[B_copy_tab[total_i_index]][j][t] / y_cond_tabsB[B_copy_tab[total_i_index]][t]);  // Weight the first and second observations together
			else
				for(j = 0; j < hmm->y_dim; j++)
					for(t = 0; t < data->seq_len; t++)
						phi[total_i_index][j][t] = gamma[total_i_index][t] *  
						(y_mu_square_exp_vecs[A_copy_tab[total_i_index]+5*B_copy_tab[total_i_index]][j][t] / 
						y_cond_tabs[A_copy_tab[total_i_index]+5*B_copy_tab[total_i_index]][t]); // Weight the first and second observations together
#else
			for(t = 0; t < data->seq_len; t++)
				phi[total_i_index][0][t] = gamma[total_i_index][t]; 
#endif
		}

	return 0;
}


// Compute the auxillary psi parameters used to re-estimate the model parameters later.
//  The meaning of Psi is: Psi[i][j][t] = Pr(X_t=i,X_{t+1}=j, Y_1,..,Y_n) (it should be conditioned on the Y's but this normalization is done elsewhere)
long ComputePsi(hmm_model *hmm, hmm_data *data, double *alpha[MAX_X_VALS], double *beta[MAX_X_VALS], 
				double *scale, double *y_cond_tabs[MAX_X_VALS], 
				double *psi[MAX_X_VALS][MAX_X_VALS])
{
	long i, j, t;

	// a simple loop ...
	if(data->miss_data == 0)   // no missing data 
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				for(t = 0; t < data->seq_len-1; t++)		
					psi[i][j][t] = hmm->M[i][j] * alpha[i][t] * beta[j][t+1] * y_cond_tabs[j][t+1];
	else
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				for(t = 0; t < data->seq_len-1; t++)		
					psi[i][j][t] = hmm->M_POWERS[data->loc_diff_vec[t+1]][i][j] * alpha[i][t] * beta[j][t+1] * y_cond_tabs[j][t+1];

	return 0; 

}


// Compute the auxillary psi parameters used to re-estimate the model parameters later, 
// in the special SNPs case
long ComputePsiSNPs(hmm_model *hmm, hmm_data *data, double *alpha[36], double *beta[36], 
				double *scale, double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS],
				double psi_sum[36][36])
{
	long i, j, t;
	long i1,i2,j1,j2,i_geno,j_geno,total_i_index,total_j_index;
	double Markov_M_Mult;

	// a simple loop ...
	if(data->miss_data == 0)   // no missing data 
		for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)	 // two genotypes indexes
			for(i1 = 0; i1 < hmm->x_dim; i1++)			// two copy numbers
				for(i2 = 0; i2 < hmm->x_dim; i2++)	
				{
					total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 
					for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	 // two genotypes indexes
						for(j1 = 0; j1 < hmm->x_dim; j1++)			// two copy numbers
							for(j2 = 0; j2 < hmm->x_dim; j2++)	
							{
								Markov_M_Mult = hmm->M[i1][j1] * hmm->M[i2][j2]; // Optimization: Do multiplication once to save time ..
								total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 						
								// New! Compute directly the psi_sums !!! 
								psi_sum[total_i_index][total_j_index] = 0; 	// Initilize sums to zero	
								if(hmm->gauss_dim == 1)
									for(t = 0; t < data->seq_len-1; t++)	
										psi_sum[total_i_index][total_j_index] += ( Markov_M_Mult * 
											hmm->place_M[BIT(i_geno,0)][BIT(j_geno,0)][t] * 
											hmm->place_M[BIT(i_geno,1)][BIT(j_geno,1)][t] *  						
											alpha[total_i_index][t] * beta[total_j_index][t+1] * 
											y_cond_tabs[A_copy_tab[total_j_index]][t+1] *  
											y_cond_tabsB[B_copy_tab[total_j_index]][t+1] );
								else
									for(t = 0; t < data->seq_len-1; t++)	
										psi_sum[total_i_index][total_j_index] += ( Markov_M_Mult * 
											hmm->place_M[BIT(i_geno,0)][BIT(j_geno,0)][t] * 
											hmm->place_M[BIT(i_geno,1)][BIT(j_geno,1)][t] *  						
											alpha[total_i_index][t] * beta[total_j_index][t+1] * 
											y_cond_tabs[A_copy_tab[total_j_index]+5*B_copy_tab[total_j_index]][t+1] ); 
							}
				}
	else // didn't fix this case yet !!!! Not working .. no ENERGY
		for(i = 0; i < 36; i++)
			for(j = 0; j < 36; j++)
				for(t = 0; t < data->seq_len-1; t++)						
					psi_sum[i][j] += hmm->M_POWERS[data->loc_diff_vec[t+1]][i][j] * alpha[i][t] * beta[j][t+1] * y_cond_tabs[j][t+1];

	return 0; 

}





// For the continuous case : Here we replace the tabs of N by computing
// Mixture of Gaussians coefficients !! 
// Here y_cond_tabs[i][t] = Pr(Y = y_vec[t] | X = i) = \sum_{j}  N[i][j] * Gaussian_{i,j} (y_vec[t]) 
// Note : If we want tommorow a different continuous distribution, we chagne only (?) here !!! 
// We compute here also for discrete distributions to gain (?) performance !!! 
long Compute_y_cond_tabs(hmm_model *hmm,  hmm_data *data, 
					double *y_cond_tabs[MAX_X_VALS], 
					double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS])
{
	long i, j, t;

	if(hmm->y_type == DISCRETE)
	{
		for(i = 0; i < hmm->x_dim; i++)
			for(t = 0; t < data->seq_len; t++)
				y_cond_tabs[i][t] = hmm->N[i][data->y_vec_int[t]];
	}
	else
	{
		if(hmm->place_flag == 0)  // standard model. No problem ...	
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					for(t = 0; t < data->seq_len; t++)
						y_mu_square_exp_vecs[i][j][t] =  ONE_OVER_2SQRT_PI * hmm->N[i][j] * hmm->SIGMA_inv[i][j] *
						exp( -(data->y_vec[t] - hmm->MU[i][j])*(data->y_vec[t] - hmm->MU[i][j]) * 0.5 *hmm->SIGMA_inv[i][j]*hmm->SIGMA_inv[i][j] ); 
		else   // here different gaussian for every place 
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					for(t = 0; t < data->seq_len; t++)
					{
						y_mu_square_exp_vecs[i][j][t] =  ONE_OVER_2SQRT_PI * hmm->N[i][j] * 
							exp( -(data->y_vec[t] - hmm->place_gaussian_mu[i][j][t])*(data->y_vec[t] - hmm->place_gaussian_mu[i][j][t]) / 
							(2.0*hmm->place_gaussian_sigma[i][j][t]*hmm->place_gaussian_sigma[i][j][t]) ) / hmm->place_gaussian_sigma[i][j][t]; 
					}
		// Now update the b-tabs according to the mixtures 
		for(i = 0; i < hmm->x_dim; i++)
		{
			// start already with the first one to save time ..
			for(t = 0; t < data->seq_len; t++)
				y_cond_tabs[i][t] = y_mu_square_exp_vecs[i][0][t];
			for(j = 1; j < hmm->y_dim; j++)
				for(t = 0; t < data->seq_len; t++)
					y_cond_tabs[i][t] += y_mu_square_exp_vecs[i][j][t];					
			for(t = 0; t < data->seq_len; t++) 	// Now make sure that they are not too small !!!!!
				y_cond_tabs[i][t] = MAX(y_cond_tabs[i][t], EPSILON);
		}
	}  // if DISCRETE

	return 0; 

}




// For the continuous case : Here we replace the tabs of N by computing
// Mixture of Gaussians coefficients !! 
// Here y_cond_tabs[i][t] = Pr(Y = y_vec[t] | X = i) = \sum_{j}  N[i][j] * Gaussian_{i,j} (y_vec[t]) 
// Note : If we want tommorow a different continuous distribution, we chagne only (?) here !!! 
// We compute here also for discrete distributions to gain (?) performance !!! 
// This is special for the SNPs model!!! 
long Compute_y_cond_tabsSNPs(hmm_model *hmm,  hmm_data *data, 
					double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS], 
					double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS], 
					double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS])
{
	long i, t;
	long B_copy, A_copy;
	double half_sigma_inv_sqr_A, half_sigma_inv_sqr_B, N_times_sigma_inv_A, N_times_sigma_inv_B;
	double tmp_x[2];
	double det_sigma; 

	if(hmm->place_flag == 0)  // standard model. No problem ...	
		for(i = 0; i < hmm->x_dim; i++)
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
			{
				half_sigma_inv_sqr_A = 0.5 * hmm->SIGMA_inv[i][j]*hmm->SIGMA_inv[i][j]; // Compute before to save time ..
				for(t = 0; t < data->seq_len; t++)
					y_mu_square_exp_vecs[i][j][t] = ONE_OVER_2SQRT_PI * hmm->N[i][j] * hmm->SIGMA_inv[i][j] *
					exp( -(data->y_vec[t] - hmm->MU[i][j])*(data->y_vec[t] - hmm->MU[i][j]) * half_sigma_inv_sqr_A); 
			}
#else
		{
			half_sigma_inv_sqr_A = 0.5 * hmm->SIGMA_inv[i][0]*hmm->SIGMA_inv[i][0]; // Compute before to save time ..
				for(t = 0; t < data->seq_len; t++)
					y_mu_square_exp_vecs[i][0][t] = ONE_OVER_2SQRT_PI * hmm->SIGMA_inv[i][0] *
					exp( -(data->y_vec[t] - hmm->MU[i][0])*(data->y_vec[t] - hmm->MU[i][0]) * half_sigma_inv_sqr_A); 
		}
#endif

	else   // here different Gaussian for every place 
	{
		if(hmm->gauss_dim == 1) // seperate one dimensional gaussians
			for(A_copy = 0; A_copy <= 4; A_copy++)
			{
#ifdef Y_DIM
				for(j = 0; j < hmm->y_dim; j++)			// the mixture chosen
				{
					B_copy = A_copy;
					half_sigma_inv_sqr_A = 0.5* hmm->SIGMA_inv[A_copy][j]*hmm->SIGMA_inv[A_copy][j];
					half_sigma_inv_sqr_B = 0.5* hmm->SIGMA_inv[B_copy][j]*hmm->SIGMA_inv[B_copy][j];
					N_times_sigma_inv_A = ONE_OVER_2SQRT_PI * hmm->N[i1][j] * hmm->SIGMA_inv[A_copy][j];
					N_times_sigma_inv_B = ONE_OVER_2SQRT_PI * hmm->N[i1][j] * hmm->SIGMA_inv[B_copy][j];						
					for(t = 0; t < data->seq_len; t++)
					{
						y_mu_square_exp_vecs[B_copy][j][t] = N_times_sigma_inv_B * 
							exp( -(data->y_vec[t] - hmm->MU[B_copy][j])*(data->y_vec[t] - hmm->MU[B_copy][j]) * 
							half_sigma_inv_sqr_B);   
						y_mu_square_exp_vecsB[A_copy][j][t] = N_times_sigma_inv_A * 
							exp( -(data->y_vecB[t] - hmm->MU[A_copy][j])*(data->y_vecB[t] - hmm->MU[A_copy][j]) * 
							half_sigma_inv_sqr_A);  
					}
				}
#else
				half_sigma_inv_sqr_A = 0.5* hmm->SIGMA_inv[A_copy][0]*hmm->SIGMA_inv[A_copy][0];
				N_times_sigma_inv_A = ONE_OVER_2SQRT_PI * hmm->SIGMA_inv[A_copy][0];
				B_copy = A_copy;
				half_sigma_inv_sqr_B = 0.5* hmm->SIGMA_inv[B_copy][0]*hmm->SIGMA_inv[B_copy][0];
				N_times_sigma_inv_B = ONE_OVER_2SQRT_PI * hmm->SIGMA_inv[B_copy][0];						
				
				for(t = 0; t < data->seq_len; t++)
				{					
					y_mu_square_exp_vecs[B_copy][0][t] = N_times_sigma_inv_B * 
							exp( -(data->y_vec[t] - hmm->MU[B_copy][0])*(data->y_vec[t] - hmm->MU[B_copy][0]) * 
							half_sigma_inv_sqr_B);   
					y_mu_square_exp_vecsB[A_copy][0][t] = N_times_sigma_inv_A * 
						exp( -(data->y_vecB[t] - hmm->MU[A_copy][0])*(data->y_vecB[t] - hmm->MU[A_copy][0]) * 
						half_sigma_inv_sqr_A);  
				}
						
#endif
			}
		else // here dimension is 2 - we also have a different gaussian for each SNP 
			for(A_copy = 0; A_copy <= 4; A_copy++)
				for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
				{
//					return 99; // early termination
// #ifdef Y_DIM currently everything is the same - we use the y dim for the 2-d mu vec and sigma matrix instead of the MoG 

					for(t = 0; t < data->seq_len; t++)
					{
						tmp_x[0] = data->y_vec[t] - hmm->place_gaussian_mu[A_copy+5*B_copy][0][t];
						tmp_x[1] = data->y_vecB[t] - hmm->place_gaussian_mu[A_copy+5*B_copy][1][t];
/***/
						tmp_x[0] = hmm->place_gaussian_sigma_inv[A_copy+5*B_copy][0][t]*tmp_x[0]*tmp_x[0] + 
							hmm->place_gaussian_sigma_inv[A_copy+5*B_copy][1][t]*tmp_x[1]*tmp_x[1] + 
							2*hmm->place_gaussian_sigma_inv[A_copy+5*B_copy][2][t]*tmp_x[0]*tmp_x[1];
/***
						tmp_x[0] = tmp_x[0]*tmp_x[0] + tmp_x[1]*tmp_x[1];

/***/
						det_sigma = hmm->place_gaussian_sigma[A_copy+5*B_copy][0][t]*hmm->place_gaussian_sigma[A_copy+5*B_copy][1][t] - 
							hmm->place_gaussian_sigma[A_copy+5*B_copy][2][t]*hmm->place_gaussian_sigma[A_copy+5*B_copy][2][t];
/***/
						y_mu_square_exp_vecs[A_copy+5*B_copy][0][t] = 
/***							1.0/(hmm->place_gaussian_sigma_inv[A_copy+5*B_copy][0][t]*tmp_x[0]*tmp_x[0]+
							hmm->place_gaussian_sigma_inv[A_copy+5*B_copy][1][t]*tmp_x[1]*tmp_x[1]); ***/ 
							exp (-tmp_x[0] / 2.0) * ONE_OVER_2PI / sqrt(det_sigma);
/****/
					}
			}



	} // end of else on gaussian dimension


//	return 0; // early termination

	// Now update the y-tabs according to the mixtures 
	if(hmm->gauss_dim == 1) // seperate one dimensional gaussians
		for(A_copy = 0; A_copy <= 4; A_copy++)		
		{
			// start already with the first one to save time ..
			for(t = 0; t < data->seq_len; t++)
			{
				y_cond_tabs[A_copy][t] = y_mu_square_exp_vecs[A_copy][0][t];
				y_cond_tabsB[A_copy][t] = y_mu_square_exp_vecsB[A_copy][0][t];
			}
#ifdef Y_DIM
			for(j = 1; j < hmm->y_dim; j++)
				for(t = 0; t < data->seq_len; t++)
				{				
					y_cond_tabs[A_copy][t] += y_mu_square_exp_vecs[A_copy][j][t];					
					y_cond_tabsB[A_copy][t] += y_mu_square_exp_vecsB[A_copy][j][t];
				}
#endif
			for(t = 0; t < data->seq_len; t++)
			{
				y_cond_tabs[A_copy][t] = MAX(y_cond_tabs[A_copy][t], 0.0001*EPSILON);
				y_cond_tabsB[A_copy][t] = MAX(y_cond_tabsB[A_copy][t], 0.0001*EPSILON);
			} 
		}
	else // here gaussian dimension is 2
		for(A_copy = 0; A_copy <= 4; A_copy++)		
			for(B_copy = 0/*MAX(1-A_copy,0)*/; (B_copy <= 4-A_copy) && (A_copy*B_copy!=3); B_copy++) // here we need both A copy and B copy together
			{
				// start already with the first one to save time ..
				for(t = 0; t < data->seq_len; t++)
				{
					y_cond_tabs[A_copy+5*B_copy][t] = y_mu_square_exp_vecs[A_copy+5*B_copy][0][t];
//					y_cond_tabsB[A_copy+5*B_copy][t] = y_mu_square_exp_vecsB[A_copy+5*B_copy][0][t];
				}
#ifdef Y_DIM
				for(j = 1; j < hmm->y_dim; j++)
					for(t = 0; t < data->seq_len; t++)
					{				
						y_cond_tabs[A_copy+5*B_copy][t] += y_mu_square_exp_vecs[A_copy+5*B_copy][j][t];					
//						y_cond_tabsB[A_copy+5*B_copy][t] += y_mu_square_exp_vecsB[A_copy+5*B_copy][j][t];
					}
#endif
				for(t = 0; t < data->seq_len; t++) // avoid too small probabilities
				{
					y_cond_tabs[A_copy+5*B_copy][t] = MAX(y_cond_tabs[A_copy+5*B_copy][t], 0.0001*EPSILON);
//					y_cond_tabsB[A_copy+5*B_copy][t] = MAX(y_cond_tabsB[A_copy+5*B_copy][t], 0.0001*EPSILON);
				} 
			}


	return 0; 

}



// Update the current model parameters to improve likelihood.
// This is the M-step of the EM algorithm.
long UpdateModelParams(hmm_model *hmm, hmm_data *data, double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS],
						 double *gamma[MAX_X_VALS], double *phi[MAX_X_VALS][MAX_Y_VALS],
						 double *psi[MAX_X_VALS][MAX_X_VALS])
{
	long i, j, t;
	double gamma_sum[MAX_X_VALS];
	double psi_sum[MAX_X_VALS][MAX_X_VALS];
	double gamma_sum_for_given_output[MAX_X_VALS][MAX_Y_VALS];
	double spare; // For correcting the M matrix

	// Compute the sums for the gammas and psi's. Note : We go only until T-1 !!
	for(i = 0; i < hmm->x_dim; i++)
	{
		// Init to zero 
		gamma_sum[i] = 0;
		for(j = 0; j < hmm->x_dim; j++)
			psi_sum[i][j] = 0; 
		// loop to sum
		for(t = 0; t < data->seq_len-1; t++)
		{
			gamma_sum[i] += gamma[i][t];
			for(j = 0; j < hmm->x_dim; j++)
				psi_sum[i][j] += psi[i][j][t];
		}
		if(hmm->y_type == DISCRETE)
		{
			for(j = 0; j < hmm->y_dim; j++)
				gamma_sum_for_given_output[i][j] = 0; 
			for(t = 0; t < data->seq_len-1; t++)
				gamma_sum_for_given_output[i][data->y_vec_int[t]] += gamma[i][t];
		}
	}
	// First update the initial distribution PI
	for(i = 0; i < hmm->x_dim; i++)
		hmm->PI[i] = gamma[i][0];
	for(i = 0; i < hmm->x_dim; i++)
		hmm->PI[i] = MAX(EPSILON, hmm->PI[i]);

	if(hmm->update_M_flag) 
	{
		// Now update the transition probability matrix M
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				hmm->M[i][j] = MAX(EPSILON, psi_sum[i][j] / gamma_sum[i]); 		// avoid too little probs !!!! 

		// avoid too large probs if given as a request !!!!
		if(hmm->use_bounds)
		{
			for(i = 0; i < hmm->x_dim; i++)
			{
				spare = 0; 
				for(j = 0; j < hmm->x_dim; j++)
				{
					spare = spare + MAX(hmm->M[i][j]-hmm->M_upperbounds[i][j], 0);
					hmm->M[i][j] = MIN(hmm->M_upperbounds[i][j], hmm->M[i][j]);
				}
				// We now need to adjust so that the rows sum up to one ! 
				for(j = 0; j < hmm->x_dim; j++)
				{
					if(spare <= (hmm->M_upperbounds[i][j]- hmm->M[i][j]))
					{
						hmm->M[i][j] += spare; 
						break;
					}
					else
					{
						spare = spare - (hmm->M_upperbounds[i][j]- hmm->M[i][j]);
						hmm->M[i][j] += (hmm->M_upperbounds[i][j]- hmm->M[i][j]);
					}
				}									
			}
		}

		// Compute powers of M if needed !!! 
		if(hmm->miss_data)
		{	
			// Copy the transition matrix 
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->x_dim; j++)
					hmm->M_POWERS[0][i][j] = hmm->M[i][j];
			// Multiply the transition matrix 
			for(t = 1; t < MAX_POWER_OF_TRANS_MATRIX; t++)
				MatrixMultiply(hmm->M, hmm->x_dim, hmm->x_dim, hmm->M_POWERS[t-1], hmm->x_dim, 
								hmm->M_POWERS[t]);
		}


	}

	// Finally update the emission matrix N
	for(i = 0; i < hmm->x_dim; i++)
		gamma_sum[i] += gamma[i][data->seq_len-1];
	if(hmm->y_type == DISCRETE)
		if(hmm->update_N_flag)
		{
			// add the last 
			for(i = 0; i < hmm->x_dim; i++)
				gamma_sum_for_given_output[i][data->y_vec_int[data->seq_len-1]] += gamma[i][data->seq_len-1];
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->N[i][j] = gamma_sum_for_given_output[i][j] / gamma_sum[i];
		}
	else  // CONTINUOUS
	{
		// Here the continuous (MoG) case : 
		double phi_sum[MAX_X_VALS][MAX_Y_VALS];
		double phi_obs_sum[MAX_X_VALS][MAX_Y_VALS];
		double phi_obs_var_sum[MAX_X_VALS][MAX_Y_VALS];

		// compute phi sum
		if(hmm->place_flag == 0) // gaussians independent of place
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
				{
					phi_sum[i][j] = 0; phi_obs_sum[i][j] = 0; phi_obs_var_sum[i][j] = 0;
					for(t = 0; t < data->seq_len-1; t++)
					{
						phi_sum[i][j] += phi[i][j][t];
						phi_obs_sum[i][j] += phi[i][j][t] * data->y_vec[t];
						phi_obs_var_sum[i][j] += phi[i][j][t] * 
							(data->y_vec[t] - hmm->MU[i][j]) * (data->y_vec[t] - hmm->MU[i][j]);
					}
				}
		else   // here everything depends on place
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
				{
					phi_sum[i][j] = 0; phi_obs_sum[i][j] = 0; phi_obs_var_sum[i][j] = 0;
					for(t = 0; t < data->seq_len-1; t++)
					{
						phi_sum[i][j] += phi[i][j][t];
						phi_obs_sum[i][j] += phi[i][j][t] * data->y_vec[t];
						phi_obs_var_sum[i][j] += phi[i][j][t] * 
							(data->y_vec[t] -  hmm->place_gaussian_mu[i][j][t]) * (data->y_vec[t] -  hmm->place_gaussian_mu[i][j][t]); 						
					}
				}

		if(hmm->update_N_flag)
		{
			// Update the N's (mixture coefficients)
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->N[i][j] = phi_sum[i][j] / gamma_sum[i];
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
		if(hmm->update_MU_flag)
		{
			// Update MU's (the means)
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->MU[i][j] = phi_obs_sum[i][j]/phi_sum[i][j];

			// make sure that the lowest level is no more than one, and the highest
			// level is no less then one. 
			if(hmm->fold_change_flag)
			{
				long max_index = 0; long min_index = 0; 
				double average_mu, max_mu, min_mu;
				max_mu = -9999999999.9;
				min_mu = 9999999999.9;

				// first find the lowest and highest level 
				for(i = 0; i < hmm->x_dim; i++)  // the X's states 
				{
					average_mu = 0;
					for(j = 0; j < hmm->y_dim; j++)
						average_mu += hmm->MU[i][j] * hmm->N[i][j];

					if(average_mu > max_mu)
					{
						max_index = i;
						max_mu = average_mu; 
					}
					if(average_mu < min_mu)
					{
						min_index = i;
						min_mu = average_mu; 
					}
				}
				// Now set to one if needed - replaced by zero !!!!!!!
				for(j = 0; j < hmm->y_dim; j++)
					hmm->MU[max_index][j] = MAX(hmm->MU[max_index][j], 0.0);   // 1.0
				for(j = 0; j < hmm->y_dim; j++)
					hmm->MU[min_index][j] = MIN(hmm->MU[min_index][j], 0.0);   // 1.0
			}
		}		
		if(hmm->update_SIGMA_flag) 			// Update SIGMA's (the standard deviations)
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->SIGMA[i][j] = MAX(EPSILON, sqrt(phi_obs_var_sum[i][j]/phi_sum[i][j])); 			// avoid too little sigma !!!!
	}  // end if continuous
	
	// avoid too little probs !!!!
	if(hmm->update_N_flag)
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
				hmm->N[i][j] = MAX(EPSILON, hmm->N[i][j]);
	ComputeModelAuxillaryParameters(hmm);  // Must be last !!!

	return 0; 
}	


// Update the current model parameters to improve likelihood.
// This is the M-step of the EM algorithm, in the SNPs version of the model
double UpdateModelParamsSNPs(hmm_model *hmm, hmm_data *data, 
						   double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS],
						   double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS],
						   double *gamma[36], double *phi[36][MAX_Y_VALS],
						   double psi_sum[36][36])
{
	long i, j, t;
	long i1,i2,j1,j2,i_geno,j_geno,total_i_index,total_j_index;

	long A_copy, B_copy; 
	double gamma_sum[36];

	double collapse_gamma_sum[MAX_X_VALS];
	double collapse_psi_sum[MAX_X_VALS][MAX_X_VALS];


	// make sure all collapse vectors are zero
	for(i=0; i<hmm->x_dim; i++)
	{
		collapse_gamma_sum[i] = 0;
		for(j=0; j<hmm->x_dim; j++)
			collapse_psi_sum[i][j] = 0;	
	}


	// Compute the sums for the gammas and psi's. Note : We go only until T-1 !!
	for(total_i_index = 0; total_i_index < 36; total_i_index++)		
	{
		// Init to zero 
		gamma_sum[total_i_index] = 0;
		for(t = 0; t < data->seq_len-1; t++)
			gamma_sum[total_i_index] += gamma[total_i_index][t];
	}
			
	// For the model parameters we need to 'collapse' the values we've calculated .. do some summation ..		
	for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)	
		for(i1 = 0; i1 < hmm->x_dim; i1++)	
			for(i2 = 0; i2 < hmm->x_dim; i2++)	
			{
				total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 	
				collapse_gamma_sum[i1] += gamma_sum[total_i_index];
				collapse_gamma_sum[i2] += gamma_sum[total_i_index];
				for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)	
					for(j1 = 0; j1 < hmm->x_dim; j1++)	
						for(j2 = 0; j2 < hmm->x_dim; j2++)	
						{
							total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 	
							collapse_psi_sum[i1][j1] += psi_sum[total_i_index][total_j_index];
							collapse_psi_sum[i2][j2] += psi_sum[total_i_index][total_j_index];
						}
			}
	for(i=0; i<hmm->x_dim; i++) // normalize to one (???) - need to check for correct normalization
	{
		collapse_gamma_sum[i] *= 0.5;
		for(j=0; j<hmm->x_dim; j++)
			collapse_psi_sum[i][j] *= 0.5;	
	}

			
	// First update the initial distribution PI
	for(i=0; i<hmm->x_dim; i++)
		hmm->PI[i] = 0;
	for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)	
		for(i1 = 0; i1 < hmm->x_dim; i1++)	
			for(i2 = 0; i2 < hmm->x_dim; i2++)	
			{
				total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 			
				hmm->PI[i1] += gamma[total_i_index][0];
				hmm->PI[i2] += gamma[total_i_index][0];
			}
	for(i=0; i<hmm->x_dim; i++)
		hmm->PI[i] = MAX(EPSILON, hmm->PI[i]/2);


	if(hmm->update_M_flag) 
	{
		// Now update the transition probability matrix M
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				hmm->M[i][j] = MAX(EPSILON, collapse_psi_sum[i][j] / collapse_gamma_sum[i]); 		// avoid too little probs !!!! 


		// avoid too large probs if given as a request !!!! (Not supprted yet!!!!!)
		if(hmm->miss_data) 		// Now compute powers of M if needed !!!  Not Working Yet !!!!! 
		{	
			// Copy the transition matrix 
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->x_dim; j++)
					hmm->M_POWERS[0][i][j] = hmm->M[i][j];
			// Multiply the transition matrix 
			for(t = 1; t < MAX_POWER_OF_TRANS_MATRIX; t++)
				MatrixMultiply(hmm->M, hmm->x_dim, hmm->x_dim, hmm->M_POWERS[t-1], hmm->x_dim, 
								hmm->M_POWERS[t]);
		}
	}

	// Finally update the emission matrix N. Start by adding the last gamma contribution
	for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)	
		for(i1 = 0; i1 < hmm->x_dim; i1++)	
			for(i2 = 0; i2 < hmm->x_dim; i2++)	
			{
				total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 			
				gamma_sum[total_i_index] += gamma[total_i_index][data->seq_len-1];
				collapse_gamma_sum[i1] += (gamma[total_i_index][data->seq_len-1]*0.5);
				collapse_gamma_sum[i2] += (gamma[total_i_index][data->seq_len-1]*0.5);
			}
	  
	// CONTINUOUS - Should Work!!! 
	// Here the continuous (MoG) case : 
	double phi_sum[36][MAX_Y_VALS];
	double phi_obs_sum[36][MAX_Y_VALS];
	double phi_obs_sumB[36][MAX_Y_VALS];
	double phi_obs_var_sum[36][MAX_Y_VALS];
	double phi_obs_var_sumB[36][MAX_Y_VALS];


	double collapse_phi_sum[MAX_X_VALS][MAX_Y_VALS];
	double collapse_phi_sumA[MAX_X_VALS][MAX_Y_VALS];
	double collapse_phi_sumB[MAX_X_VALS][MAX_Y_VALS];
	double collapse_phi_obs_sum[MAX_X_VALS][MAX_Y_VALS];
	double collapse_phi_obs_sumB[MAX_X_VALS][MAX_Y_VALS];
	double collapse_phi_obs_var_sum[MAX_X_VALS][MAX_Y_VALS];
	double collapse_phi_obs_var_sumB[MAX_X_VALS][MAX_Y_VALS];
	

	// compute phi sum
	if( (hmm->place_flag == 0) || (hmm->place_flag == 1)) // gaussians independent of place, eventhough flag is one (for place_M)
		for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)	
			for(i1 = 0; i1 < hmm->x_dim; i1++)	
				for(i2 = 0; i2 < hmm->x_dim; i2++)	
				{
					total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 			
					A_copy = (1^BIT(i_geno,0)) * i1 + (1^BIT(i_geno,1)) * i2; // A copy number
					B_copy = BIT(i_geno,0) * i1 + BIT(i_geno,1) * i2;	// B copy number
#ifdef Y_DIM
					for(j = 0; j < hmm->y_dim; j++)
					{						
						phi_sum[total_i_index][j] = 0; 
						phi_obs_sum[total_i_index][j] = 0; phi_obs_sumB[total_i_index][j] = 0; 
						phi_obs_var_sum[total_i_index][j] = 0; phi_obs_var_sumB[total_i_index][j] = 0;
						for(t = 0; t < data->seq_len-1; t++)
						{
							phi_sum[total_i_index][j] += phi[total_i_index][j][t];
							phi_obs_sum[total_i_index][j] += phi[total_i_index][j][t] * data->y_vec[t];
							phi_obs_sumB[total_i_index][j] += phi[total_i_index][j][t] * data->y_vecB[t];
							phi_obs_var_sum[total_i_index][j] += phi[total_i_index][j][t] * 
								(data->y_vec[t] - hmm->MU[A_copy][j]) * (data->y_vec[t] - hmm->MU[A_copy][j]);
							phi_obs_var_sumB[total_i_index][j] += phi[total_i_index][j][t] * 
								(data->y_vecB[t] - hmm->MU[B_copy][j]) * (data->y_vecB[t] - hmm->MU[B_copy][j]);								
						}
					}
#else
					phi_sum[total_i_index][0] = 0; 
					phi_obs_sum[total_i_index][0] = 0; phi_obs_sumB[total_i_index][0] = 0; 
					phi_obs_var_sum[total_i_index][0] = 0; phi_obs_var_sumB[total_i_index][0] = 0;
					for(t = 0; t < data->seq_len-1; t++)
					{
						phi_sum[total_i_index][0] += phi[total_i_index][0][t];
						phi_obs_sum[total_i_index][0] += phi[total_i_index][0][t] * data->y_vec[t];
						phi_obs_sumB[total_i_index][0] += phi[total_i_index][0][t] * data->y_vecB[t];
						phi_obs_var_sum[total_i_index][0] += phi[total_i_index][0][t] * 
							(data->y_vec[t] - hmm->MU[A_copy][0]) * (data->y_vec[t] - hmm->MU[A_copy][0]);
						phi_obs_var_sumB[total_i_index][0] += phi[total_i_index][0][t] * 
							(data->y_vecB[t] - hmm->MU[B_copy][0]) * (data->y_vecB[t] - hmm->MU[B_copy][0]);								
					}
#endif
				}
	else   // here everything depends on place
		for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)	
			for(i1 = 0; i1 < hmm->x_dim; i1++)	
				for(i2 = 0; i2 < hmm->x_dim; i2++)	
				{
					total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 			
#ifdef Y_DIM					
					for(j = 0; j < hmm->y_dim; j++)
					{
						phi_sum[total_i_index][j] = 0; phi_obs_sum[total_i_index][j] = 0; phi_obs_sumB[total_i_index][j] = 0; phi_obs_var_sum[total_i_index][j] = 0;
						for(t = 0; t < data->seq_len-1; t++)
						{
							phi_sum[total_i_index][j] += phi[total_i_index][j][t];
							phi_obs_sum[total_i_index][j] += phi[total_i_index][j][t] * data->y_vec[t];
							phi_obs_var_sum[total_i_index][j] += phi[total_i_index][j][t] * 
								(data->y_vec[t] -  hmm->place_gaussian_mu[total_i_index][j][t]) * 
								(data->y_vec[t] -  hmm->place_gaussian_mu[total_i_index][j][t]); 						
						}
					}
#else
					phi_sum[total_i_index][0] = 0; phi_obs_sum[total_i_index][0] = 0; phi_obs_sumB[total_i_index][0] = 0; phi_obs_var_sum[total_i_index][0] = 0;
					for(t = 0; t < data->seq_len-1; t++)
					{
						phi_sum[total_i_index][0] += phi[total_i_index][0][t];
						phi_obs_sum[total_i_index][0] += phi[total_i_index][0][t] * data->y_vec[t];
						phi_obs_var_sum[total_i_index][0] += phi[total_i_index][0][t] * 
							(data->y_vec[t] -  hmm->place_gaussian_mu[total_i_index][0][t]) * 
							(data->y_vec[t] -  hmm->place_gaussian_mu[total_i_index][0][t]); 						
					}
#endif
				}


	// deal with the 'collapsed' sums: 
	for(i = 0; i < (2*hmm->x_dim-1); i++)			
#ifdef Y_DIM
		for(j = 0; j < hmm->y_dim; j++)
		{
			collapse_phi_sum[i][j] = 0;
			collapse_phi_sumA[i][j] = 0;
			collapse_phi_sumB[i][j] = 0;
			collapse_phi_obs_sum[i][j] = 0;
			collapse_phi_obs_sumB[i][j] = 0;
			collapse_phi_obs_var_sum[i][j] = 0;
			collapse_phi_obs_var_sumB[i][j] = 0;
		}
#else
	{
		collapse_phi_sum[i][0] = 0;
		collapse_phi_sumA[i][0] = 0;
		collapse_phi_sumB[i][0] = 0;
		collapse_phi_obs_sum[i][0] = 0;
		collapse_phi_obs_sumB[i][0] = 0;
		collapse_phi_obs_var_sum[i][0] = 0;
		collapse_phi_obs_var_sumB[i][0] = 0;
	}
#endif

	for(i_geno = 0; i_geno < (hmm->x_dim2*hmm->x_dim2); i_geno++)	
		for(i1 = 0; i1 < hmm->x_dim; i1++)	
			for(i2 = 0; i2 < hmm->x_dim; i2++)	
			{
				total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 					
				A_copy = (1^BIT(i_geno,0)) * i1 + (1^BIT(i_geno,1)) * i2; // A copy number
				B_copy = BIT(i_geno,0) * i1 + BIT(i_geno,1) * i2;	// B copy numnber
#ifdef Y_DIM
				for(j = 0; j < hmm->y_dim; j++)
				{
					collapse_phi_sum[i1][j] += (phi_sum[total_i_index][j]*0.5);
					collapse_phi_sum[i2][j] += (phi_sum[total_i_index][j]*0.5); // why only these are divided by two??? 
					collapse_phi_sumA[A_copy][j] += phi_sum[total_i_index][j];
					collapse_phi_sumB[B_copy][j] += phi_sum[total_i_index][j];
					collapse_phi_obs_sum[A_copy][j] += phi_obs_sum[total_i_index][j];
					collapse_phi_obs_sumB[B_copy][j] += phi_obs_sumB[total_i_index][j];
					collapse_phi_obs_var_sum[A_copy][j] += phi_obs_var_sum[total_i_index][j];
					collapse_phi_obs_var_sumB[B_copy][j] += phi_obs_var_sumB[total_i_index][j];
				}
#else
				collapse_phi_sum[i1][0] += (phi_sum[total_i_index][0]*0.5);
				collapse_phi_sum[i2][0] += (phi_sum[total_i_index][0]*0.5); // why only these are divided by two??? 
				collapse_phi_sumA[A_copy][0] += phi_sum[total_i_index][0];
				collapse_phi_sumB[B_copy][0] += phi_sum[total_i_index][0];
				collapse_phi_obs_sum[A_copy][0] += phi_obs_sum[total_i_index][0];
				collapse_phi_obs_sumB[B_copy][0] += phi_obs_sumB[total_i_index][0];
				collapse_phi_obs_var_sum[A_copy][0] += phi_obs_var_sum[total_i_index][0];
				collapse_phi_obs_var_sumB[B_copy][0] += phi_obs_var_sumB[total_i_index][0];
#endif
			}

	
	double norm_gamma_one = 0.0;
	for(i=0; i<(2*hmm->x_dim-1); i++)
#ifdef Y_DIM
		for(j=0; j<hmm->y_dim; j++)	
			norm_gamma_one += collapse_phi_sumB[i][j];
#else
		norm_gamma_one += collapse_phi_sumB[i][0];
#endif

			
	if(hmm->update_N_flag)
	{
		// Update the N's (mixture coefficients)
		for(i = 0; i < (2*hmm->x_dim-1); i++)	
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
				hmm->N[i][j] = MAX(EPSILON,  collapse_phi_sum[i][j] / (collapse_gamma_sum[i])); 	// avoid too little probs !!!!
#else
			hmm->N[i][0] = MAX(EPSILON,  collapse_phi_sum[i][0] / (collapse_gamma_sum[i])); 	// avoid too little probs !!!!
#endif

		// Normalize by sum (should not be neccessary, just extra-security)
#ifdef Y_DIM
		double sumsum;
		for(i = 0; i < (2*hmm->x_dim-1); i++)	
		{
			sumsum = 0; 
			for(j = 0; j < hmm->y_dim; j++)
				sumsum += hmm->N[i][j];
			for(j = 0; j < hmm->y_dim; j++)
				hmm->N[i][j] /= sumsum;
		}
#endif
	}
	if(hmm->update_MU_flag)
	{
		// Update MU's (the means)
		for(i = 0; i < (2*hmm->x_dim-1); i++)	
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
				hmm->MU[i][j] = (collapse_phi_obs_sum[i][j] + collapse_phi_obs_sumB[i][j]) /
				(collapse_phi_sumA[i][j] + collapse_phi_sumB[i][j]);
#else
			hmm->MU[i][0] = (collapse_phi_obs_sum[i][0] + collapse_phi_obs_sumB[i][0]) /
				(collapse_phi_sumA[i][0] + collapse_phi_sumB[i][0]);
#endif

	}		
	if(hmm->update_SIGMA_flag)
	{
		// Update SIGMA's (the standard deviations)
		for(i = 0; i < (2*hmm->x_dim-1); i++)	   // the X's states 
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
			{
				hmm->SIGMA[i][j] = sqrt( (collapse_phi_obs_var_sum[i][j] + collapse_phi_obs_var_sumB[i][j]) / 
					(collapse_phi_sumA[i][j] + collapse_phi_sumB[i][j]));  
				hmm->SIGMA[i][j] = MAX(MIN_SIGMA/*EPSILON*/, hmm->SIGMA[i][j]); // avoid too little sigma !!!!
			}
#else
		{
			hmm->SIGMA[i][0] = sqrt( (collapse_phi_obs_var_sum[i][0] + collapse_phi_obs_var_sumB[i][0]) / 
				(collapse_phi_sumA[i][0] + collapse_phi_sumB[i][0]));  
			hmm->SIGMA[i][0] = MAX(MIN_SIGMA, hmm->SIGMA[i][0]); // avoid too little sigma !!!!
		}

#endif
	}
	// end continuous


	
	// Here we need to sort them in a 'special way' by the MU's 
	PermuteModelIncreasingMean(hmm, data);  // permute EVERY time !!! 


	// We also want to keep the MU's in the desired range
	for(i = 0; i < (2*hmm->x_dim-1); i++)
		hmm->MU[i][0] = MIN( MAX(hmm->MU[i][0], MIN_MU_VAL[i]), MAX_MU_VAL[i]);
	

	ComputeModelAuxillaryParameters(hmm);  // Must be last !!!


	return norm_gamma_one; // 0; 
			
}	








// Get a smart starting point for the hmm parameters ..
// Not implemented yet ..
long SmartStartParams(hmm_model *hmm, double *y_vec, long seq_len)
{
	return 0;
}



// perform the EM algorithm to establish the model parameters from data.
// The input variales are:
// hmm - contains many flags saying what kind of model are we using (discrete/continuous observations,
//	are their place dependent transition matrices etc.)
// data - the observational data we learn from.
// max_iters - number of EM iterations to perform
// num_starts - number of different starting points for the EM algorithm
// tolerance - the amount of improvement we require in each EM iteration - otherwise we stop.
//
// The output is stored in the structure hmm, and also the log-likelihood score 
// is given in best_model_score.
long TrainModelEM(hmm_model *hmm, hmm_data *data, long max_iters, long num_starts, double tolerance, 
				  double *best_model_score)
{
	
	printf("START EM INSIDE ! miss date is %ld \n", hmm->miss_data);
	
	long i, j, k;


	// First Initilize model, use the same function so we have here a 'self-reference paradox' ..
	long x_dim = hmm->x_dim;
	long y_dim = hmm->y_dim;
	long y_type = hmm->y_type;
	long place_flag = hmm->place_flag;

	hmm->cum_flag = 0;
	hmm->log_flag = 1;


	// Now start the iterative EM ;-)
	double *y_cond_tabs[MAX_X_VALS];
	double *alpha[MAX_X_VALS];
	double *beta[MAX_X_VALS];
	double *gamma[MAX_X_VALS];
	double *phi[MAX_X_VALS][MAX_Y_VALS];
	double *psi[MAX_X_VALS][MAX_X_VALS];	// The estimated transition probs.
	double *scale = new double[data->seq_len];   // scaling to avoid underflows
	double *scale_cum = new double[data->seq_len];   // cumulative sum of the scale
	double *scale_exp = new double[data->seq_len];   // exp of the scale (ratio)
	double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS];    // (y-mu)^2


	long best_model_index = -1;


	// First allocate memory 
	for(i = 0; i < hmm->x_dim; i++)
	{
		y_cond_tabs[i] = new double[data->seq_len];
		alpha[i] = new double[data->seq_len];
		beta[i] = new double[data->seq_len];
		gamma[i] = new double[data->seq_len];
	}
	for(i = 0; i < hmm->x_dim; i++)
		for(j = 0; j < hmm->x_dim; j++)
			psi[i][j] = new double[data->seq_len];
	if(hmm->y_type == CONTINUOUS)
	{
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->y_dim; j++)
			{
				phi[i][j] = new double[data->seq_len];	
				y_mu_square_exp_vecs[i][j] = new double[data->seq_len];
			}
	}



	*best_model_score = -99999999999;  // the best score to keep ! 
	hmm_model best_hmm[MAX_STARTING_POINTS];   // here many models !! 	


	double cur_scores[MAX_STARTING_POINTS];
	double prev_scores[MAX_STARTING_POINTS];

	long iter;

	for(i = 0; i < num_starts; i++)  // try different starting points 
	{
		// Currently we use random parameters start !!!
		InitilizeModel(&(best_hmm[i]), x_dim, y_dim, y_type,  hmm->use_bounds, place_flag, hmm->special_models_flag, hmm->miss_data, hmm->fold_change_flag); 

		// Copy the upperbounds
		if(hmm->use_bounds)
			for(k = 0; k < (&(best_hmm[i]))->x_dim; k++)		
				for(j = 0; j < (&(best_hmm[i]))->x_dim; j++)
					(&(best_hmm[i]))->M_upperbounds[k][j] = hmm->M_upperbounds[k][j];

		// Copy everything
		if(place_flag)
		{
					// allocate memory ... 
			for(k = 0; k < (&(best_hmm[i]))->x_dim; k++)		
				for(j = 0; j < (&(best_hmm[i]))->y_dim; j++)
				{
					(&(best_hmm[i]))->place_gaussian_mu[k][j] = new double[data->seq_len];
					(&(best_hmm[i]))->place_gaussian_sigma[k][j] = new double[data->seq_len];
				}

			CopyModel(hmm, &(best_hmm[i]), data);
		}

		// Currently we use random parameters start !!!
		InitilizeModel(&(best_hmm[i]), x_dim, y_dim, y_type, hmm->use_bounds, place_flag, 
			hmm->special_models_flag, hmm->miss_data, hmm->fold_change_flag); 	
		// Copy the upperbounds again
		if(hmm->use_bounds)
			for(k = 0; k < (&(best_hmm[i]))->x_dim; k++)		
				for(j = 0; j < (&(best_hmm[i]))->x_dim; j++)
					(&(best_hmm[i]))->M_upperbounds[k][j] = hmm->M_upperbounds[k][j];

		// Here maybe use a more clever starting point ... 
		iter = 0;		
		cur_scores[i] = -99999999999; 
		prev_scores[i] = -99999999999.9 - 2*tolerance;
		
		// Do first pass and see if we increase ..
		while(  (iter < FIRST_ITERS) &&  (ABS(cur_scores[i] - prev_scores[i])) > tolerance  )
		{
		
			prev_scores[i] = cur_scores[i];  // update the previous
				
			// E-step : estimate probabilities	
// #define CHECK_TIMING
#ifdef CHECK_TIMING
			clock_t EM_start_time = clock();
#endif
			Compute_y_cond_tabs(&(best_hmm[i]), data, y_cond_tabs, y_mu_square_exp_vecs); 
#ifdef CHECK_TIMING
			clock_t EM_end_time = clock();
			printf("Compute_y_cond_tabs (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			forward(&(best_hmm[i]), data, y_cond_tabs, 
				alpha, scale, scale_cum, scale_exp, &(cur_scores[i]));
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Forward (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			backward(&(best_hmm[i]), data, scale_exp, y_cond_tabs, 
				beta);		
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Backward (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			ComputeMarginalGamma(&(best_hmm[i]), data, y_mu_square_exp_vecs, 
				alpha, beta, scale_exp, y_cond_tabs, 
				gamma, phi);
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Marginal (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			ComputePsi(&(best_hmm[i]), data, alpha, beta, scale, y_cond_tabs, psi);
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Psi (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			// M-step : update model parameters
			UpdateModelParams(&(best_hmm[i]), data, y_mu_square_exp_vecs, gamma, phi, psi);	
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Update (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif

			if(iter%10 == 9)
				printf("FINISHED ITER %ld SCORE : %lf  ----\n", iter, cur_scores[i]); // PrintModel(&best_hmm); printf("\n---\n Now Realy Finished Iter \n----\n");
			iter++;

		}

		printf("Did Phase 1 Start %ld Itrs %ld score %lf prev %lf best %lf\n", 
			i, iter, cur_scores[i], prev_scores[i], *best_model_score);



		if(  (cur_scores[i] > (*best_model_score)) || // here we have improved !! Run again to improve iterations !!!! 
			( ABS( (cur_scores[i] - (*best_model_score)) / (*best_model_score) ) < IMPROVEMENT_FRACTION ) ) // Not improved but close !!  Run again to improve iterations !!!!
		{
			while(  (iter < max_iters) && (ABS(cur_scores[i] - prev_scores[i]) > tolerance)  )
			{
			
				prev_scores[i] = cur_scores[i];  // update the previous
				
				// E-step : estimate probabilities	
	// #define CHECK_TIMING
#ifdef CHECK_TIMING
				clock_t EM_start_time = clock();
#endif
				Compute_y_cond_tabs(&(best_hmm[i]), data, y_cond_tabs, y_mu_square_exp_vecs); 
#ifdef CHECK_TIMING
				clock_t EM_end_time = clock();
				printf("Compute_y_cond_tabs (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				forward(&(best_hmm[i]), data, y_cond_tabs, alpha, scale, scale_cum, scale_exp, &(cur_scores[i]));
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Forward (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				backward(&(best_hmm[i]), data, scale_exp, y_cond_tabs, beta);		
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Backward (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				ComputeMarginalGamma(&(best_hmm[i]), data, y_mu_square_exp_vecs, alpha, beta, scale_exp, y_cond_tabs, 
					gamma, phi);
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Marginal (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				ComputePsi(&(best_hmm[i]), data, alpha, beta, scale, y_cond_tabs, psi);
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Psi (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				// M-step : update model parameters
				UpdateModelParams(&(best_hmm[i]), data, y_mu_square_exp_vecs, gamma, phi, psi);	
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Update (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif

				if(iter%10 == 9)
					printf("FINISHED ITER %ld SCORE : %lf  ----\n", iter, cur_scores[i]); // PrintModel(&best_hmm); printf("\n---\n Now Realy Finished Iter \n----\n");
				iter++;

			}

			printf("Did Phase 2 Start %ld Itrs %ld score %lf prev %lf\n", i, iter, cur_scores[i], prev_scores[i]);

			// Now we need to condition since we can accept also below best score
			if(cur_scores[i] > (*best_model_score)) 
			{
				*best_model_score = cur_scores[i];
				best_model_index = i;  // save the place from which came the best model !! 
			}
		}   // end if of phase 2


	}   // end loop on starting points 


//////////////////////////////////////////////////////////////////////////////
	// Now do Phase 3 : Take the best offer from all starting points and improve it. Do another max_iters iterations 

	iter = max_iters;

	printf("BEFORE STAGE 3 : best_score %lf prev_score %lf diff %lf\n", 
		cur_scores[best_model_index], prev_scores[best_model_index], 
		cur_scores[best_model_index] - prev_scores[best_model_index]);
	while(  (iter < 2*max_iters) && (ABS(cur_scores[best_model_index] - prev_scores[best_model_index]) > tolerance)  )
	{			
		prev_scores[best_model_index] = cur_scores[best_model_index];  // update the previous
			
	
		// E-step : estimate probabilities	
	// #define CHECK_TIMING
#ifdef CHECK_TIMING
		clock_t EM_start_time = clock();
#endif
		Compute_y_cond_tabs(&(best_hmm[best_model_index]), data, y_cond_tabs, y_mu_square_exp_vecs); 
#ifdef CHECK_TIMING
		clock_t EM_end_time = clock();
		printf("Compute_y_cond_tabs (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
		forward(&(best_hmm[best_model_index]), data, y_cond_tabs, alpha, scale, scale_cum, scale_exp, &(cur_scores[best_model_index]));
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Forward (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
	
		
		backward(&(best_hmm[best_model_index]), data, scale_exp, y_cond_tabs, beta);		
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Backward (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
		ComputeMarginalGamma(&(best_hmm[best_model_index]), data, y_mu_square_exp_vecs, alpha, beta, scale_exp, y_cond_tabs, 
				gamma, phi);
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Marginal (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
		ComputePsi(&(best_hmm[best_model_index]), data, alpha, beta, scale, y_cond_tabs, psi);
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Psi (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
		// M-step : update model parameters
		UpdateModelParams(&(best_hmm[best_model_index]), data, y_mu_square_exp_vecs, 
			gamma, phi, psi);	
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Update (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif

		if(iter%10 == 9)
			printf("FINISHED ITER %ld SCORE : %lf  ----\n", iter, cur_scores[best_model_index]); // PrintModel(&best_hmm); printf("\n---\n Now Realy Finished Iter \n----\n");
		iter++;
	}
	
	printf("Did Final Phase 3 Start %ld Itrs %ld score %lf prev %lf\n", 
		best_model_index, iter, cur_scores[best_model_index], prev_scores[best_model_index]);
			
		// Now we need to condition since we can accept also below best score
	*best_model_score = cur_scores[best_model_index];
	CopyModel(&(best_hmm[best_model_index]), hmm, data);

	printf("Inside EM : Best MODEL IS : \n");
	PrintModel(&(best_hmm[best_model_index]));
	printf("-----------\n\n");


////////////////////////////////////////////////////////////////////////////



	// return all the cums !!! 
	hmm->cum_flag = 1; PermuteModelIncreasingMean(hmm, data); 

	Compute_y_cond_tabs(hmm, data, y_cond_tabs,
		y_mu_square_exp_vecs); 
	forward(hmm, data, y_cond_tabs, alpha, scale, scale_cum, scale_exp, 
		best_model_score);


	// free all allocated buffers
	delete scale;
	delete scale_cum;
	delete scale_exp;
	for(i = 0; i < hmm->x_dim; i++)
	{
		delete y_cond_tabs[i];
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
				delete phi[i][j];
				delete y_mu_square_exp_vecs[i][j];
			}
	}
	printf("Finished EM. Best Score is %lf\n", *best_model_score);
	return 0;
}








// perform the EM algorithm to establish the model parameters from data,
// on the special HMM SNPs model ! 
double TrainModelEMSNPs(hmm_model *hmm, hmm_data *data, long max_iters, long num_starts, double tolerance, 
				  double *best_model_score, double *tmp_vec)
{	
//	printf("START EM SNPs INSIDE ! miss data is %ld \n", hmm->miss_data);
	
	long i, j, k;
	long  total_i_index; 

	// First Initilize model, use the same function so we have here a 'self-reference paradox' ..
	long x_dim = hmm->x_dim;
	long y_dim = hmm->y_dim;
	long y_type = hmm->y_type;
	long place_flag = hmm->place_flag;

	hmm->cum_flag = 0;
	hmm->log_flag = 1;
	hmm->seq_len = data->seq_len;

	// Now start the iterative EM ;-)
	double *y_cond_tabs[MAX_2D_X_VALS];
	double *y_cond_tabsB[MAX_2D_X_VALS];
	double *alpha[36];
	double *beta[36];
	double *gamma[36];
	double *phi[36][MAX_Y_VALS];
//	double *psi[36][36];	// The estimated transition probs.
	double psi_sum[36][36]; // The sums of the estimated transition probs.
	double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS];    // (y-mu)^2
	double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS];    // (y-mu)^2


//////////////////////////////////////////////////////////////////////////////  Phase 0 ////////////////
/****************************************************************************/

	double *scale = new double[data->seq_len];   // scaling to avoid underflows
	double *scale_cum = new double[data->seq_len];   // cumulative sum of the scale
	double *scale_exp = new double[data->seq_len];   // exp of the scale (ratio)
/****************************************************************************/
//////////////////////////////////////////////////////////////////////////////  Ended Phase 0 ////////////////
	long best_model_index = -1;

	// First allocate memory 
	long max_x_states = 36;
	*best_model_score = -99999999999;  // the best score to keep ! 
	hmm_model best_hmm[MAX_STARTING_POINTS];   // here many models !! 	
	double cur_scores[MAX_STARTING_POINTS];
	double prev_scores[MAX_STARTING_POINTS];

	long iter;



//////////////////////////////////////////////////////////////////////////////  Phase 0 ////////////////
/****************************************************************************/
	for(total_i_index = 0; total_i_index < 36; total_i_index++)		
	{
		// Seperate according to whether the source or destination are zero ..
		alpha[total_i_index] = new double[data->seq_len];
		beta[total_i_index] = new double[data->seq_len];
		gamma[total_i_index] = new double[data->seq_len];
		if(total_i_index <= 4)
		{
			y_cond_tabs[total_i_index] = new double[data->seq_len];
			y_cond_tabsB[total_i_index] = new double[data->seq_len];
		}

/***
		for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)		// state at time t-1
			for(j1 = 0; j1 < hmm->x_dim; j1++)		// state at time t-1
				for(j2 = 0; j2 < hmm->x_dim; j2++)		// state at time t-1
						{
							// Seperate according to whether the source or destination are zero ..
							total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
							psi[total_i_index][total_j_index] = new double[data->seq_len]; // seems that we allocated it fine here
						}
***/
		if(hmm->y_type == CONTINUOUS)
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
			{
				phi[total_i_index][j] = new double[data->seq_len];	
				if(total_i_index <= 4)
				{
					y_mu_square_exp_vecs[total_i_index][j] = new double[data->seq_len];
					y_mu_square_exp_vecsB[total_i_index][j] = new double[data->seq_len];
				}
			}
#else
		{
			phi[total_i_index][0] = new double[data->seq_len];	
			if(total_i_index <= 4)
			{
				y_mu_square_exp_vecs[total_i_index][0] = new double[data->seq_len];
				y_mu_square_exp_vecsB[total_i_index][0] = new double[data->seq_len];
			}
		}
#endif
	}
/****************************************************************************/
//////////////////////////////////////////////////////////////////////////////  Ended Phase 0 ////////////////

//	tmp_vec[1] = hmm->x_dim;


	for(i = 0; i < num_starts; i++)  // try different starting points 
	{
		// Currently we use random parameters start !!!

		InitilizeModel(&(best_hmm[i]), x_dim, y_dim, y_type,  hmm->use_bounds, place_flag, 
			hmm->special_models_flag, hmm->miss_data, hmm->fold_change_flag);  
		(&(best_hmm[i]))->update_M_flag = 1; (&(best_hmm[i]))->update_N_flag = 0; // we can save a little time by not touching N 
		(&(best_hmm[i]))->update_MU_flag =1; (&(best_hmm[i]))->update_SIGMA_flag = 1;   // Try not to update and see if we improve

		// Copy the upperbounds
		if(hmm->use_bounds)
			for(k = 0; k < (&(best_hmm[i]))->x_dim; k++)		
				for(j = 0; j < (&(best_hmm[i]))->x_dim; j++)
					(&(best_hmm[i]))->M_upperbounds[k][j] = hmm->M_upperbounds[k][j];
		// Copy everything
		if(place_flag)
		{
			// allocate memory ... 
			for(k = 0; k < 2; k++)		
				for(j = 0; j < 2; j++)	
				{
					(&(best_hmm[i]))->place_M[k][j] = new double[data->seq_len];
					(&(best_hmm[i]))->place_M_cum[k][j] = new double[data->seq_len];
				}
			CopyModel(hmm, &(best_hmm[i]), data);
		}

		// Currently we use random parameters start !!!
		InitilizeModel(&(best_hmm[i]), x_dim, y_dim, y_type, hmm->use_bounds, place_flag, 
			hmm->special_models_flag, hmm->miss_data, hmm->fold_change_flag); 	
		(&(best_hmm[i]))->update_M_flag = 1; (&(best_hmm[i]))->update_N_flag = 0; // we can save a little time by not touching N
		(&(best_hmm[i]))->update_MU_flag = 1; (&(best_hmm[i]))->update_SIGMA_flag = 1; // Try not to update and see if we improve
		// Copy the upperbounds
		if(hmm->use_bounds)
			for(k = 0; k < (&(best_hmm[i]))->x_dim; k++)		
				for(j = 0; j < (&(best_hmm[i]))->x_dim; j++)
					(&(best_hmm[i]))->M_upperbounds[k][j] = hmm->M_upperbounds[k][j];



		// Temporary: Try giving a 'close to optimal' starting point
//		best_hmm[i].MU[0][0] = 1.04;  best_hmm[i].MU[1][0] = 2.003; best_hmm[i].MU[2][0] = 1.99; best_hmm[i].MU[3][0] = 4.01; best_hmm[i].MU[4][0] = 4.9; 
//		best_hmm[i].SIGMA[0][0] = 0.1;  best_hmm[i].SIGMA[1][0] = 0.17; best_hmm[i].SIGMA[2][0] = 0.03; best_hmm[i].SIGMA[3][0] = 0.07; best_hmm[i].SIGMA[4][0] = 0.13; 


		// Here maybe use a more clever starting point ... 
		iter = 0;		
		cur_scores[i] = -99999999999; 
		prev_scores[i] = -99999999999.9 - 2*tolerance; // make gap between them large enough
		
//////////////////////////////////////////////////////////////////////////////  Phase 1 ////////////////
/****************************************************************************/



		// Do first pass and see if we increase ..
		while(  (iter < FIRST_ITERS) && (cur_scores[i] - prev_scores[i]) > tolerance  )
		{
		
			prev_scores[i] = cur_scores[i];  // update the previous
			// E-step : estimate probabilities	
// #define CHECK_TIMING
#ifdef CHECK_TIMING
			clock_t EM_start_time = clock();
#endif

			Compute_y_cond_tabsSNPs(&(best_hmm[i]), data, y_cond_tabs, y_cond_tabsB,
				 y_mu_square_exp_vecs, y_mu_square_exp_vecsB); // last repeat (dummy)

#ifdef CHECK_TIMING
			clock_t EM_end_time = clock();
			printf("Compute_y_cond_tabs (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			forwardSNPs(&(best_hmm[i]), data, y_cond_tabs, y_cond_tabsB,  
				alpha, scale, scale_cum, scale_exp, &(cur_scores[i]));
//			tmp_vec[i] = cur_scores[i]; // save the scores always
			best_model_index = i; // try to finish the loop 			

#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Forward SNPs (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
//			return 1977; // still 'no fly'
			backwardSNPs(&(best_hmm[i]), data, scale_exp, y_cond_tabs, y_cond_tabsB,
				beta);		

//////////////////////////////////////////////////////////////////////////////  Part of Phase 1 ////////////////
/****************************************************************************/


#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Backward SNPs (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			ComputeMarginalGammaSNPs(&(best_hmm[i]), data, y_mu_square_exp_vecs, y_mu_square_exp_vecsB, 
				alpha, beta, scale_exp, y_cond_tabs, y_cond_tabsB,
				gamma, phi);

#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Marginal (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			ComputePsiSNPs(&(best_hmm[i]), data, alpha, beta, scale, y_cond_tabs,  y_cond_tabsB, psi_sum);
#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Psi (sec.) : %lf\n", 
				double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			// M-step : update model parameters				
			UpdateModelParamsSNPs(&(best_hmm[i]), data, y_mu_square_exp_vecs, y_mu_square_exp_vecsB, 
				gamma, phi, psi_sum);	
/***
			if(iter <= 1)
			{
				for(j=0; j < 5; j++)
					tmp_vec[50+20*iter+j] = best_hmm[i].SIGMA_inv[j][0]; 
				for(j=0; j < 5; j++)
					tmp_vec[50+20*iter+j+5] = best_hmm[i].SIGMA[j][0]; 
				for(j=0; j < 5; j++)
					tmp_vec[50+20*iter+j+10] = best_hmm[i].MU[j][0]; 

				// Print the collapse phi values : (are they zeros???) 
				for(j=0; j<5; j++)
				{
					tmp_vec[50+20*iter+j+15] = 100;
					for(k=0; k<36; k++)
						tmp_vec[50+20*iter+j+15] += y_cond_tabs[k][j]; // gamma[k][j]; //  phi[k][0][j];
				}
			}
***/


#ifdef CHECK_TIMING
			EM_end_time = clock();
			printf("Update (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
			EM_start_time = clock();
#endif
			if(iter%10 == 9)
				printf("FINISHED ITER %ld SCORE : %lf  ----\n", iter, cur_scores[i]); // PrintModel(&best_hmm); printf("\n---\n Now Realy Finished Iter \n----\n");
			iter++;


/****************************************************************************/
//////////////////////////////////////////////////////////////////////////////  Ended Part of Phase 1 ////////////////


		}

		printf("Did Phase 1 Start %ld Itrs %ld score %lf prev %lf best %lf\n", 
			i, iter, cur_scores[i], prev_scores[i], *best_model_score);
/****************************************************************************/
//////////////////////////////////////////////////////////////////////////////  Ended Phase 1 ////////////////


//////////////////////////////////////////////////////////////////////////////  Phase 2 ////////////////
/****************************************************************************/
		if(  (cur_scores[i] > (*best_model_score)) || // here we have improved !! Run again to improve iterations !!!! 
			( ABS( (cur_scores[i] - (*best_model_score)) / (*best_model_score) ) < IMPROVEMENT_FRACTION ) ) // Not improved but close enough !!  Run again to improve iterations !!!!
		{

			while(  (iter < max_iters) && ((cur_scores[i] - prev_scores[i]) > tolerance)  )
			{
			
				prev_scores[i] = cur_scores[i];  // update the previous
				
				// E-step : estimate probabilities	
	// #define CHECK_TIMING
#ifdef CHECK_TIMING
				clock_t EM_start_time = clock();
#endif
				Compute_y_cond_tabsSNPs(&(best_hmm[i]), data, y_cond_tabs, y_cond_tabsB, y_mu_square_exp_vecs, y_mu_square_exp_vecsB); // last repeat (dummy)
#ifdef CHECK_TIMING
				clock_t EM_end_time = clock();
				printf("Compute_y_cond_tabs (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				forwardSNPs(&(best_hmm[i]), data, y_cond_tabs, y_cond_tabsB, alpha, scale, scale_cum, scale_exp, &(cur_scores[i]));
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Forward SNPs (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				backwardSNPs(&(best_hmm[i]), data, scale_exp, y_cond_tabs, y_cond_tabsB, beta);		
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Backward (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				ComputeMarginalGammaSNPs(&(best_hmm[i]), data, y_mu_square_exp_vecs, y_mu_square_exp_vecsB, 
					alpha, beta, scale_exp, y_cond_tabs,  y_cond_tabsB, 
					gamma, phi);
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Marginal (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				ComputePsiSNPs(&(best_hmm[i]), data, alpha, beta, scale, y_cond_tabs, y_cond_tabsB, psi_sum);
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Psi (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				// M-step : update model parameters
				UpdateModelParamsSNPs(&(best_hmm[i]), data, y_mu_square_exp_vecs,  y_mu_square_exp_vecsB, gamma, phi, psi_sum);	
#ifdef CHECK_TIMING
				EM_end_time = clock();
				printf("Update (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
				EM_start_time = clock();
#endif
				if(iter%10 == 9)
					printf("FINISHED ITER %ld SCORE : %lf  ----\n", iter, cur_scores[i]); // PrintModel(&best_hmm); printf("\n---\n Now Realy Finished Iter \n----\n");
				iter++;
			}

			printf("Did Phase 2 Start %ld Itrs %ld score %lf prev %lf\n", i, iter, cur_scores[i], prev_scores[i]);

			// Now we need to condition since we can accept also below best score
			if(cur_scores[i] > *best_model_score) 
			{
				*best_model_score = cur_scores[i];
				best_model_index = i;  // save the place from which came the best model !! 
			}
		}   // end if of phase 2


/****************************************************************************/
//////////////////////////////////////////////////////////////////////////////

	}   // end loop on starting points 


	// save best model and score
//	tmp_vec[8] = best_model_index; 
//	tmp_vec[9] = (*best_model_score);

//////////////////////////////////////////////////////////////////////////////
/****************************************************************************/
	// Now do Phase 3 : Take the best offer from all starting points and improve it. Do another max_iters iterations 
	iter = max_iters;
	

	printf("BEFORE STAGE 3 : best_score %lf prev_score %lf diff %lf\n", 
		cur_scores[best_model_index], prev_scores[best_model_index], 
		cur_scores[best_model_index] - prev_scores[best_model_index]);
	while(  (iter < 2*max_iters) && ((cur_scores[best_model_index] - prev_scores[best_model_index]) > tolerance)  )
	{			
		prev_scores[best_model_index] = cur_scores[best_model_index];  // update the previous
			
		// E-step : estimate probabilities	
	// #define CHECK_TIMING
#ifdef CHECK_TIMING
		clock_t EM_start_time = clock();
#endif
		Compute_y_cond_tabsSNPs(&(best_hmm[best_model_index]), data, y_cond_tabs, y_cond_tabsB, y_mu_square_exp_vecs, y_mu_square_exp_vecsB); // last repeat (dummy)
#ifdef CHECK_TIMING
		clock_t EM_end_time = clock();
		printf("Compute_y_cond_tabs (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
		forwardSNPs(&(best_hmm[best_model_index]), data, y_cond_tabs, y_cond_tabsB, alpha, scale, scale_cum, scale_exp, &(cur_scores[best_model_index]));
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Forward SNPs (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
	
		
		backwardSNPs(&(best_hmm[best_model_index]), data, scale_exp, y_cond_tabs, y_cond_tabsB, beta);		
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Backward SNPs (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
		ComputeMarginalGammaSNPs(&(best_hmm[best_model_index]), data, 
			y_mu_square_exp_vecs, y_mu_square_exp_vecsB, alpha, beta, scale_exp, y_cond_tabs, y_cond_tabsB,
				gamma, phi);
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Marginal (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
		ComputePsiSNPs(&(best_hmm[best_model_index]), data, alpha, beta, scale, y_cond_tabs, y_cond_tabsB, psi_sum);
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Psi (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif
		// M-step : update model parameters
		UpdateModelParamsSNPs(&(best_hmm[best_model_index]), data, y_mu_square_exp_vecs, y_mu_square_exp_vecsB, 
			gamma, phi, psi_sum);	
#ifdef CHECK_TIMING
		EM_end_time = clock();
		printf("Update (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
		EM_start_time = clock();
#endif

		if(iter%10 == 9)
			printf("FINISHED ITER %ld SCORE : %lf  ----\n", iter, cur_scores[best_model_index]); // PrintModel(&best_hmm); printf("\n---\n Now Realy Finished Iter \n----\n");
		iter++;
	}
	


	printf("Did Final Phase 3 Start %ld Itrs %ld score %lf prev %lf\n", 
		best_model_index, iter, cur_scores[best_model_index], prev_scores[best_model_index]);
			
		// Now we need to condition since we can accept also below best score
	*best_model_score = cur_scores[best_model_index];


/****************************************************************************/
////////////////////////////////////////////////////////////////////////////


	CopyModel(&(best_hmm[best_model_index]), hmm, data);
//	tmp_vec[7] = best_model_index;
	
	//	return best_model_index;


	printf("Inside EM : Best MODEL IS : \n");
	PrintModel(&(best_hmm[best_model_index]));
	printf("-----------\n\n");



	// return all the cums !!! 
	hmm->cum_flag = 1; ComputeModelAuxillaryParameters(hmm); ///////PermuteModelIncreasingMean(hmm, data); 

////////////////////////////////////////////////////////////////////////////
/****************************************************************************/

	Compute_y_cond_tabsSNPs(hmm, data, y_cond_tabs, y_cond_tabsB,
		y_mu_square_exp_vecs, y_mu_square_exp_vecsB);  // last repeat (dummy)
	forwardSNPs(hmm, data, y_cond_tabs, y_cond_tabsB,  alpha, scale, scale_cum, scale_exp, 
		best_model_score);

/****************************************************************************/
////////////////////////////////////////////////////////////////////////////


//	printf("At the end of all EM AFTER permuting the HMM SCORE IS %lf\n", *best_model_score);
//	PrintModel(hmm);

  
	// free all allocated buffers
	delete scale;
	delete scale_cum;
	delete scale_exp;


	for(total_i_index = 0; total_i_index < 36; total_i_index++)		
	{
		// Seperate according to whether the source or destination are zero ..
		delete alpha[total_i_index];
		delete beta[total_i_index];
		delete gamma[total_i_index];
		if(total_i_index <= 4)
		{
			delete y_cond_tabs[total_i_index];
			delete y_cond_tabsB[total_i_index];
		}

/****
				for(j_geno = 0; j_geno < (hmm->x_dim2*hmm->x_dim2); j_geno++)		// state at time t-1
					for(j1 = 0; j1 < hmm->x_dim; j1++)		// state at time t-1
						for(j2 = 0; j2 < hmm->x_dim; j2++)		// state at time t-1
						{
							// Seperate according to whether the source or destination are zero ..
							total_j_index = multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]; 
							delete psi[total_i_index][total_j_index]; // seems that we allocate it fine here
						}
***/
		if(hmm->y_type == CONTINUOUS)
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
			{
				delete phi[total_i_index][j];	
				if(total_i_index <= 4)
				{
					delete y_mu_square_exp_vecs[total_i_index][j];
					delete y_mu_square_exp_vecsB[total_i_index][j];
				}
			}
#else
		{
			delete phi[total_i_index][0];	
			if(total_i_index <= 4)
			{
				delete y_mu_square_exp_vecs[total_i_index][0];
				delete y_mu_square_exp_vecsB[total_i_index][0];
			}
		}
#endif
	}

	for(i = 0; i < num_starts; i++)  // delete models for the different starting points 
		if(place_flag)
		{
			// delete memory ... 
			for(k = 0; k < 2; k++)		
				for(j = 0; j < 2; j++)		
				{
					delete (&(best_hmm[i]))->place_M[k][j];
					delete (&(best_hmm[i]))->place_M_cum[k][j];
				}
		}

  /****************************************************************************/
////////////////////////////////////////////////////////////////////////////

	printf("Finished EM SNPs. Best Score is %lf\n", *best_model_score);
	return tmp_vec[5]; /// hmm->seq_len;
}







// Copy all models parameters. We copy from hmm1 to hmm2
long CopyModel(hmm_model *hmm1, hmm_model *hmm2, hmm_data *data)
{
	long i, j, t, k;

	// copy types 
	hmm2->x_dim = hmm1->x_dim;
	hmm2->x_dim2 = hmm1->x_dim2;
	hmm2->y_dim = hmm1->y_dim;
	hmm2->y_type = hmm1->y_type;
	hmm2->place_flag = hmm1->place_flag;
	hmm2->use_bounds = hmm1->use_bounds;
	hmm2->special_models_flag = hmm1->special_models_flag;
//	hmm2->seq_len = hmm1->seq_len;

	
	// copy parameters
	for(i = 0; i < hmm2->x_dim; i++)
		hmm2->PI[i] = hmm1->PI[i];	
	for(i = 0; i < hmm2->x_dim; i++)
		for(j = 0; j < hmm2->x_dim; j++)
			hmm2->M[i][j] = hmm1->M[i][j];
	if(hmm2->use_bounds)
		for(i = 0; i < hmm2->x_dim; i++)
			for(j = 0; j < hmm2->x_dim; j++)
				hmm2->M_upperbounds[i][j] = hmm1->M_upperbounds[i][j];

	if(hmm1->special_models_flag)
	{
		for(i = 0; i < 2*hmm2->x_dim-1; i++)
#ifdef Y_DIM
			for(j = 0; j < hmm2->y_dim; j++)
				hmm2->N[i][j] = hmm1->N[i][j];
#else
			hmm2->N[i][0] = hmm1->N[i][0];
#endif
		if(hmm2->y_type == CONTINUOUS)
			for(i = 0; i < 2*hmm2->x_dim-1; i++)
#ifdef Y_DIM
				for(j = 0; j < hmm2->y_dim; j++)
				{
					hmm2->MU[i][j] = hmm1->MU[i][j];
					hmm2->SIGMA[i][j] = hmm1->SIGMA[i][j];
				}
#else
			{
				hmm2->MU[i][0] = hmm1->MU[i][0];
				hmm2->SIGMA[i][0] = hmm1->SIGMA[i][0];
			}
#endif
	}
	else
	{
		for(i = 0; i < hmm2->x_dim; i++)
			for(j = 0; j < hmm2->y_dim; j++)
				hmm2->N[i][j] = hmm1->N[i][j];
		if(hmm2->y_type == CONTINUOUS)
			for(i = 0; i < hmm2->x_dim; i++)
				for(j = 0; j < hmm2->y_dim; j++)
				{
					hmm2->MU[i][j] = hmm1->MU[i][j];
					hmm2->SIGMA[i][j] = hmm1->SIGMA[i][j];
				}
	}


	// copy the means and sigmas
	if(hmm1->place_flag)
	{
		if(hmm1->special_models_flag == 0)
			for(i = 0; i < hmm2->x_dim; i++)
				for(j = 0; j < hmm2->y_dim; j++)
					for(t = 0; t < data->seq_len; t++)
					{
						hmm2->place_gaussian_mu[i][j][t] = hmm1->place_gaussian_mu[i][j][t];
						hmm2->place_gaussian_sigma[i][j][t] = hmm1->place_gaussian_sigma[i][j][t];
					}
		else
			for(i = 0; i < hmm1->x_dim2; i++)		
				for(j = 0; j < hmm2->x_dim2; j++)	
					for(t = 0; t < data->seq_len; t++)
						hmm2->place_M[i][j][t] = hmm1->place_M[i][j][t]; 
	}

	if(data->miss_data)  // copy the power matrices
	{
		for(k = 0; k < MAX_POWER_OF_TRANS_MATRIX; k++)
		for(i = 0; i < hmm2->x_dim; i++)
			for(j = 0; j < hmm2->x_dim; j++)
				hmm2->M_POWERS[k][i][j] = hmm1->M_POWERS[k][i][j];

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

	printf("X Transition Matrix :\n");
	for(i = 0; i < hmm->x_dim; i++)
		PrintDoubleVec(hmm->M[i], hmm->x_dim);

	printf("X CUMCUM Transition Matrix :\n");
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


// permute the indexes of the model, thus enabling a more fair 
// comparison of the trained and the original model !!!
long PermuteModelIncreasingMean(hmm_model *hmm, hmm_data *data)
{
	long i, j, t, k; 
	double temp, summy;
	double mean_vec[MAX_X_VALS+MAX_Y_VALS];	
	long indexes[MAX_X_VALS];
	double helperx[MAX_X_VALS][MAX_X_VALS];
	double helpery[MAX_X_VALS][MAX_Y_VALS];
	long max_x_states = hmm->x_dim;
	if(hmm->special_models_flag)
		max_x_states = 2*hmm->x_dim-1;

	if(hmm->y_type == DISCRETE)  // here permute according to conditional means, assuming the Y's are 0,1,..
	{
		for(i = 0; i < hmm->x_dim; i++)
		{
			mean_vec[i] = 0; 
			for(j = 0; j < hmm->y_dim; j++)
				mean_vec[i] += hmm->N[i][j] * j; 
		}
	}
	else	// here permute according to the conditional means
	{
		if(  (hmm->place_flag == 0)  || hmm->special_models_flag )// not according to place 
			for(i = 0; i < max_x_states; i++)
			{

#ifdef Y_DIM
				mean_vec[i] = 0; 
				for(j = 0; j < hmm->y_dim; j++)
					mean_vec[i] += hmm->N[i][j] * hmm->MU[i][j];
#else
				mean_vec[i] = hmm->MU[i][0];
#endif
			}
		else					// here permute according to place of average 
			for(i = 0; i < hmm->x_dim; i++)
			{
				mean_vec[i] = 0; 
				for(j = 0; j < hmm->y_dim; j++)
					for(t = 0; t < data->seq_len; t++)
						mean_vec[i] += hmm->N[i][j] * hmm->place_gaussian_mu[i][j][t];
				mean_vec[i] /= data->seq_len; // normalize to a reasonable size 
			}
	}

	



	if(hmm->special_models_flag == 0)
	{
		// Now sort according to mean !! 	
		DoQuicksort(mean_vec, max_x_states, indexes);


		// Change the initial distribution
		for(i = 0; i < hmm->x_dim; i++)
			helperx[i][0] = hmm->PI[indexes[i]];
		for(i = 0; i < hmm->x_dim; i++)
			hmm->PI[i] = helperx[i][0];
		// change the X's M's variables 
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				helperx[i][j] = hmm->M[indexes[i]][indexes[j]];
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				hmm->M[i][j] = helperx[i][j];
		// Permute the powers matrices 
		if(data->miss_data)
			for(k = 0; k < MAX_POWER_OF_TRANS_MATRIX; k++)
			{
				for(i = 0; i < hmm->x_dim; i++)
					for(j = 0; j < hmm->x_dim; j++)
						helperx[i][j] = hmm->M_POWERS[k][indexes[i]][indexes[j]];
			
				for(i = 0; i < hmm->x_dim; i++)
					for(j = 0; j < hmm->x_dim; j++)
						hmm->M_POWERS[k][i][j] = helperx[i][j];
			}
	}
	else // here the special flag is on
	{
		for(i = 0; i < hmm->x_dim; i++) // here we must make M 'diagonal-like'
		{
			DoQuicksort(hmm->M[i], hmm->x_dim, indexes);
			for(j = 0; j < hmm->x_dim; j++)
				helperx[i][j] = hmm->M[i][j];
			hmm->M[i][i] = helperx[i][hmm->x_dim-1]; k=0;
			for(j = 0; j < hmm->x_dim; j++)
				if(j != i) 
					hmm->M[i][j] = helperx[i][k++];
		}
	
		
		// Here make the diagonal elements big
		for(i = 0; i < hmm->x_dim; i++) // here we must make M 'diagonal-like'
		{
			hmm->M[i][i] = MAX(hmm->M[i][i],0.9); summy = 0.0;
			for(j = 0; j < hmm->x_dim; j++)
				if(j != i)
					summy += hmm->M[i][j];
			for(j = 0; j < hmm->x_dim; j++)
				if(j != i)
					hmm->M[i][j] *= ((1-hmm->M[i][i]) / summy);
		}



		double pi_vec[MAX_X_VALS];

		// Here we must make the second column the greatest
		MatrixStationaryVec(hmm->M, hmm->x_dim, pi_vec);
		DoQuicksort(pi_vec, hmm->x_dim, indexes);

		if(indexes[hmm->x_dim-1] != 1) // we want the biggest to be second. Swap two rows and columns
		{
			for(i = 0; i < hmm->x_dim; i++) 
			{
				temp = hmm->M[i][indexes[hmm->x_dim-1]];
				hmm->M[i][indexes[hmm->x_dim-1]] = hmm->M[i][1];
				hmm->M[i][1] = temp;
			}
			for(i = 0; i < hmm->x_dim; i++) 
			{
				temp = hmm->M[indexes[hmm->x_dim-1]][i];
				hmm->M[indexes[hmm->x_dim-1]][i] = hmm->M[1][i];
				hmm->M[1][i] = temp;
			}
		}



		

		// We also want to force the stationary distribution to spent most of its time on the 'middle' level.
		temp = hmm->M[0][1]; hmm->M[0][1] = MAX(temp, hmm->M[0][2]); hmm->M[0][2] = MIN(temp, hmm->M[0][2]);
		temp = hmm->M[2][1]; hmm->M[2][1] = MAX(temp, hmm->M[2][0]); hmm->M[2][0] = MIN(temp, hmm->M[2][0]);



		// Now sort according to mean !! 	
		DoQuicksort(mean_vec, max_x_states, indexes);

	}


	// Change the Y's N's variables 
	for(i = 0; i < max_x_states; i++)
#ifdef Y_DIM
		for(j = 0; j < hmm->y_dim; j++)
			helpery[i][j] = hmm->N[indexes[i]][j];
#else
		helpery[i][0] = hmm->N[indexes[i]][0];
#endif
	for(i = 0; i < max_x_states; i++)
#ifdef Y_DIM
		for(j = 0; j < hmm->y_dim; j++)
			hmm->N[i][j] = helpery[i][j];
#else
		hmm->N[i][0] = helpery[i][0];
#endif


	if(hmm->y_type == CONTINUOUS)  // Change MU and SIGMA
	{
		for(i = 0; i < max_x_states; i++)
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
				helpery[i][j] = hmm->MU[indexes[i]][j];
#else
			helpery[i][0] = hmm->MU[indexes[i]][0];
#endif
		for(i = 0; i < max_x_states; i++)
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
				hmm->MU[i][j] = helpery[i][j];
#else
			hmm->MU[i][0] = helpery[i][0];
#endif
		for(i = 0; i < max_x_states; i++)
#ifdef Y_DIM
			for(j = 0; j < hmm->y_dim; j++)
				helpery[i][j] = hmm->SIGMA[indexes[i]][j];
#else
			helpery[i][0] = hmm->SIGMA[indexes[i]][0];
#endif

			for(i = 0; i < max_x_states; i++)
#ifdef Y_DIM
				for(j = 0; j < hmm->y_dim; j++)
					hmm->SIGMA[i][j] = helpery[i][j];
#else
				hmm->SIGMA[i][0] = helpery[i][0];
#endif


		// Change the MU and SIGMA of places 
		if(hmm->special_models_flag == 0)
		{
			if(hmm->place_flag)
			{
				printf("SWAP PLACES !!!\n");
				for(t = 0; t < data->seq_len; t++)  // go over all of them and change 
				{
					for(i = 0; i < hmm->x_dim; i++)
#ifdef Y_DIM
						for(j = 0; j < hmm->y_dim; j++)
							helpery[i][j] = hmm->place_gaussian_mu[indexes[i]][j][t];
#else
						helpery[i][0] = hmm->place_gaussian_mu[indexes[i]][0][t];
#endif
						
					for(i = 0; i < hmm->x_dim; i++)
#ifdef Y_DIM
						for(j = 0; j < hmm->y_dim; j++)
							hmm->place_gaussian_mu[indexes[i]][j][t] = helpery[i][j];
#else
						hmm->place_gaussian_mu[indexes[i]][0][t] = helpery[i][0];
#endif


					for(i = 0; i < hmm->x_dim; i++)
#ifdef Y_DIM
						for(j = 0; j < hmm->y_dim; j++)
							helpery[i][j] = hmm->place_gaussian_sigma[indexes[i]][j][t];
#else
						helpery[i][0] = hmm->place_gaussian_sigma[indexes[i]][0][t];
#endif
	
					for(i = 0; i < hmm->x_dim; i++)
#ifdef Y_DIM
						for(j = 0; j < hmm->y_dim; j++)
							hmm->place_gaussian_sigma[indexes[i]][j][t] = helpery[i][j];
#else
						hmm->place_gaussian_sigma[indexes[i]][0][t] = helpery[i][0];
#endif

				}
			}
			// Now permute also the mixtures components !! Since they are hidden, they might 
			// not be ordered in the same order !!! 
			else
			{
				for(i = 0; i < hmm->x_dim; i++)
				{ 
					// copy the MU components 
					for(j = 0; j < hmm->y_dim; j++)
						mean_vec[j] = hmm->MU[i][j];

					printf("mean_mix vec : %lf %lf\n", mean_vec[0], mean_vec[1]);
				
					// Perform the sorting 
					DoQuicksort(mean_vec, hmm->y_dim, indexes);

					printf("after sort mean_mix vec : %lf %lf\n", mean_vec[0], mean_vec[1]);
					printf("indexes : %ld %ld\n", indexes[0], indexes[1]);

					// apply sort on N, MU and SIGMA
					for(j = 0; j < hmm->y_dim; j++)
						helpery[i][j] = hmm->N[i][indexes[j]];
					for(j = 0; j < hmm->y_dim; j++)
						hmm->N[i][j] = helpery[i][j];

					for(j = 0; j < hmm->y_dim; j++)
						helpery[i][j] = hmm->MU[i][indexes[j]];
					for(j = 0; j < hmm->y_dim; j++)
						hmm->MU[i][j] = helpery[i][j];

					for(j = 0; j < hmm->y_dim; j++)
						helpery[i][j] = hmm->SIGMA[i][indexes[j]];
					for(j = 0; j < hmm->y_dim; j++)
						hmm->SIGMA[i][j] = helpery[i][j];
				}
			}
		}
	}

//		return 0;

	// Change everything else 
	ComputeModelAuxillaryParameters(hmm);

	return 0;

}



// Initilize the data struct with data from outside 
long InitilizeData( hmm_data *data, long seq_len, long y_type, long miss_data, double dont_miss_prob, long copy_data, 
				   double *y_vec, long *y_vec_int, long *loc_vec)
{
	long  t;


	data->seq_len = seq_len; 
	data->miss_data = miss_data;
	data->y_type = y_type; 
	data->dont_miss_prob = dont_miss_prob;



	// Copy the needed data 
	if(copy_data)
	{
		if(data->y_type == DISCRETE)
			for(t = 0; t < data->seq_len; t++)
				data->y_vec_int[t] = y_vec_int[t];
		else
			for(t = 0; t < data->seq_len; t++)
				data->y_vec[t] = y_vec[t];
		
		// copy the locations 
		if(data->miss_data)
		{
			for(t = 0; t < data->seq_len; t++)
				data->loc_vec[t] = loc_vec[t];
			data->loc_diff_vec[0] = 0;
			for(t = 1; t < data->seq_len; t++)
				data->loc_diff_vec[t] = loc_vec[t]-loc_vec[t-1]-1;
		}
	}		

	printf("Finish Init Data \n");

	return 0; 

}



// Initilize model with given sizes. Unless stated otherwise, parameters are randomized
long InitilizeModel(hmm_model *hmm, long x_dim, long y_dim, long y_type, long use_bounds, long place_flag, 
					long special_models_flag, long miss_data, long fold_change_flag)
{
	long i, j, t;
	double sumsum, spare;
	long indexes[MAX_X_VALS];
	double mean_vec[MAX_X_VALS];
	double helpery[MAX_X_VALS][MAX_Y_VALS];



	// Copy parameters
	hmm->x_dim = x_dim;
	hmm->y_dim = y_dim; 
	hmm->y_type = y_type; 
	hmm->place_flag = place_flag;
	hmm->special_models_flag = special_models_flag;
	//	data->seq_len = seq_len; 
	hmm->miss_data = miss_data;
	hmm->fold_change_flag = fold_change_flag;
	hmm->use_bounds = use_bounds;

	hmm->update_M_flag = 1; 
	hmm->update_N_flag = 1;
	
	
	if(hmm->place_flag == 0)
	{
		hmm->x_dim2 = 2; 
		hmm->update_MU_flag = 1; 
		hmm->update_SIGMA_flag = 1;
	}
	else  // We keep the mu's fixed and do not learn them from the data
	{
		if(hmm->special_models_flag) // what is the sequence length?
		{
			hmm->x_dim2 = 2; 
			hmm->update_MU_flag = 1; 
			hmm->update_SIGMA_flag = 1;

			/*** No need to update the place_M matrix!!! ***
			for(t = 0; t < hmm->seq_len-1; t++)
				for(i = 0; i < 2; i++)
				{
					for(j = 0; j < 2; j++)
						hmm->place_M[i][j][t] = randr();
					sumsum = hmm->place_M[i][0][t] + hmm->place_M[i][1][t];
					for(j = 0; j < 2; j++)
						hmm->place_M[i][j][t] = hmm->place_M[i][j][t] / sumsum;
				}		
			***/
		}
		else
		{
			hmm->update_MU_flag = 0; 
			hmm->update_SIGMA_flag = 0;
		}
	}


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

		// Now set constraint
		if(hmm->use_bounds)
		{
			for(i = 0; i < hmm->x_dim; i++)
			{
				spare = 0; 
				for(j = 0; j < hmm->x_dim; j++)
				{
					spare = spare + MAX(hmm->M[i][j]-hmm->M_upperbounds[i][j], 0);
					hmm->M[i][j] = MIN(hmm->M_upperbounds[i][j], hmm->M[i][j]);
				}

				// We now need to adjust so that the rows sum up to one ! 
				for(j = 0; j < hmm->x_dim; j++)
				{
					if(spare <= (hmm->M_upperbounds[i][j]- hmm->M[i][j]))
					{
						hmm->M[i][j] += spare; 
						break;
					}
					else
					{
						spare = spare - (hmm->M_upperbounds[i][j]- hmm->M[i][j]);
						hmm->M[i][j] += (hmm->M_upperbounds[i][j]- hmm->M[i][j]);
					}
				}									

			}

		}
	}

	long small_noise = 0; 
	if(small_noise)
	{
		double noise = 0.05;

			// Emission matrices
		for(i = 0; i < x_dim; i++)
		{
			hmm->N[i][i] = 1-noise; // a big prob. 
			sumsum = 0; 
			for(j = 0; j < y_dim; j++)
				if(j != i)
				{
					hmm->N[i][j] = randr();
					sumsum += hmm->N[i][j];
				}
			for(j = 0; j < y_dim; j++)
				if(j!=i)
					hmm->N[i][j] = hmm->N[i][j]*noise/sumsum;  // Normalize to make sum equal one !
		}
 	}
	else
	{
		// Emission matrices
		if(hmm->special_models_flag)
			for(i = 0; i < (2*x_dim-1); i++)
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
		else
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
	}

	// MU and SIGMA. Currently  MU ~ U[-1,1] and SIGMA ~ U[0,1]
	if(hmm->y_type == CONTINUOUS)
	{
		if(hmm->special_models_flag)
		{
			for(i = 0; i < (2*x_dim-1); i++)
			{
				mean_vec[i] = 0; 
				for(j = 0; j < y_dim; j++)
				{
					hmm->MU[i][j] = 1*(2*randr()-1+1); 
					hmm->SIGMA[i][j] = randr()*1; // start with wide SIGMAs  
					mean_vec[i] += hmm->N[i][j] * hmm->MU[i][j];
				}

			}
			// Sort MU's in increasing order
			DoQuicksort(mean_vec, 2*x_dim-1, indexes);

			for(i = 0; i < (2*x_dim-1); i++)
				for(j = 0; j < hmm->y_dim; j++)
					helpery[i][j] = hmm->MU[indexes[i]][j];
			for(i = 0; i < (2*x_dim-1); i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->MU[i][j] = helpery[i][j];


		}
		else
			for(i = 0; i < x_dim; i++)
				for(j = 0; j < y_dim; j++)
				{
					hmm->MU[i][j] = 1*(2*randr()-1+1); 
					hmm->SIGMA[i][j] = randr()/2; 
				}
	}


	// Now compute powers of M if needed !!! 
	if(hmm->miss_data)
	{	
		// Copy the transition matrix 
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				hmm->M_POWERS[0][i][j] = hmm->M[i][j];


		// Multiply the transition matrix 
		for(t = 1; t < MAX_POWER_OF_TRANS_MATRIX; t++)
			MatrixMultiply(hmm->M, hmm->x_dim, hmm->x_dim, hmm->M_POWERS[t-1], hmm->x_dim, 
							hmm->M_POWERS[t]);

	}


	ComputeModelAuxillaryParameters(hmm);


	return 0;
}


// here we only use the model parameters and compute from them the helpfull cumsum and log tables
long ComputeModelAuxillaryParameters(hmm_model *hmm)
{

	long i, j, k, p;


	long max_x_states = hmm->x_dim;

	// Decide which mode are we in
	if(hmm->special_models_flag == SNPS_ONE_SAMPLE) // here work on the SNP chips
		max_x_states = 2*hmm->x_dim-1; // change the maximal value we loop to in the SNPs case		

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
		}
		if(hmm->miss_data)  // update the cumsum of powers of M
			for(k = 0; k < MAX_POWER_OF_TRANS_MATRIX; k++)		
			{
				for(i = 0; i < hmm->x_dim; i++)
				{	
					hmm->M_POWERS_cum[k][i][0] = hmm->M_POWERS[k][i][0];
					for(j = 1; j < hmm->x_dim; j++)
						hmm->M_POWERS_cum[k][i][j] = hmm->M_POWERS_cum[k][i][j-1] + 
											hmm->M_POWERS[k][i][j];
				}
			}	

/***********/
		if(hmm->special_models_flag == SNPS_ONE_SAMPLE)
		{
			if(hmm->place_flag) // here assume that we are in the SNPs case, the dim is 2
				for(p = 0; p < hmm->seq_len-1; p++)
					for(i = 0; i < hmm->x_dim2; i++)
					{
						hmm->place_M_cum[i][0][p] = hmm->place_M[i][0][p];
						for(j = 1; j < hmm->x_dim2; j++)
							hmm->place_M_cum[i][j][p] = hmm->place_M_cum[i][j-1][p] + hmm->place_M[i][j][p]; 
					}

			for(i = 0; i < (2*hmm->x_dim-1); i++)
			{
				hmm->N_cum[i][0] = hmm->N[i][0];
				for(j = 1; j < hmm->y_dim; j++)
					hmm->N_cum[i][j] = hmm->N_cum[i][j-1] + hmm->N[i][j];
			}
		}
		else
			for(i = 0; i < hmm->x_dim; i++)
			{
				hmm->N_cum[i][0] = hmm->N[i][0];
				for(j = 1; j < hmm->y_dim; j++)
					hmm->N_cum[i][j] = hmm->N_cum[i][j-1] + hmm->N[i][j];
			}
/**********/
  }


	if(hmm->log_flag)
	{
		// Compute logs for efficiency
		for(i = 0; i < hmm->x_dim; i++)
				hmm->PI_log[i] = log(hmm->PI[i]);
		for(i = 0; i < hmm->x_dim; i++)
			for(j = 0; j < hmm->x_dim; j++)
				hmm->M_log[i][j] = log(hmm->M[i][j]);
		for(i = 0; i < max_x_states; i++)
			for(j = 0; j < hmm->y_dim; j++)
				hmm->N_log[i][j] = log(hmm->N[i][j]);
		if(hmm->y_type == CONTINUOUS)
			for(i = 0; i < max_x_states; i++)
				for(j = 0; j < hmm->y_dim; j++)
					hmm->SIGMA_log[i][j] = log(hmm->SIGMA[i][j]);
	}
	if(hmm->y_type == CONTINUOUS)
		for(i = 0; i < max_x_states; i++)
			for(j = 0; j < hmm->y_dim; j++)
				hmm->SIGMA_inv[i][j] = 1.0 / hmm->SIGMA[i][j];

	return 0; 
}



// Give various measures of (dis)similarity between two HMMs
long CompareTwoModels(hmm_model *hmm1, hmm_model *hmm2)
{
	long i, j;


	double M_ssq = 0; double N_ssq = 0; 
	double MU_ssq = 0, SIGMA_ssq = 0;   // only for continuous
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
					MU_ssq += (hmm1->MU[i][j] - hmm2->MU[i][j])*(hmm1->MU[i][j] - hmm2->MU[i][j]);   // The means 
					SIGMA_ssq += (hmm1->SIGMA[i][j] - hmm2->SIGMA[i][j])*(hmm1->SIGMA[i][j] - hmm2->SIGMA[i][j]);   // The variances
					N_ssq += (hmm1->N[i][j] - hmm2->N[i][j])*(hmm1->N[i][j] - hmm2->N[i][j]);   // The mixture coefficients
				}
			}
		}
	}


	// Now do the Kullblack-Leibler Comparison
	double KL_dist1_to_2, KL_dist2_to_1;
	long iters = 20;
	hmm_data data;

	data.y_type = hmm1->y_type;
	data.seq_len = 10000;

	ComputeKLDistance(hmm1, hmm2, &data, iters, &KL_dist1_to_2);
	ComputeKLDistance(hmm2, hmm1, &data, iters, &KL_dist2_to_1);


	// Summury of comparison

	printf("Two Models Comparison : \nM_SSQ %lf   N_SSQ %lf\n", M_ssq, N_ssq);
	if(hmm1->y_type  == CONTINUOUS)
		printf("MU_SSQ %lf   SIGMA_SSQ %lf\n", MU_ssq, SIGMA_ssq);
	printf("KL12 : %lf KL21 : %lf KL_SIM : %lf\n", KL_dist1_to_2, KL_dist2_to_1, (KL_dist1_to_2+KL_dist2_to_1)/2);


	return 0; // everything is OK

}



// Simulate a sequence from a given model. We assume that the vector is already allocated .
// We output also the x vector if someone wants it ...
long SimulateSequenceFromModel(hmm_model *hmm, hmm_data *data)
{

	long j, t; // t denotes time .....

	long mix_choose;
	long B_copy, A_copy;
	long interval;

	// first generate the x'th ...
	// Choose the 1st state according to the initial distribution 
	double randy = randr();	
	if(hmm->special_models_flag == SNPS_ONE_SAMPLE)
	{
		data->x_vec[0] = 0; // see if this makes things better, i.e. never BIG numbers
		for(j = 0; j < hmm->x_dim; j++)
			if(randy <= hmm->PI_cum[j])
			{
				data->x_vec[0] = (j << 1);
				break;
			}		
		randy = randr();	
		for(j = 0; j < hmm->x_dim; j++)
			if(randy <= hmm->PI_cum[j])
			{
				data->x_vec[0] ^= (j << 5);
				break;
			}		
	}
	else
	{
		for(j = 0; j < hmm->x_dim; j++)
			if(randy <= hmm->PI_cum[j])
			{
				data->x_vec[0] = j;
				break;
			}
	}
		
	if(data->miss_data)
	{
		data->loc_vec[0]=1; // first place is one ...
		data->loc_diff_vec[0]=0; // no gap before first place ...
	}


	// Choose all other states. assume geometric distribution of intervals (that is, each place has
	// equal prob. to be missing)
	if(data->miss_data == 0)  // no missing data
	{
		if(hmm->place_flag == 1) 
		{
			if(hmm->special_models_flag == SNPS_ONE_SAMPLE)
			{
				for(t = 1; t < data->seq_len; t++)
				{
					// determine the genotype (binary data, stored in the lsb) for alpha
					randy = randr();
					for(j = 0; j < hmm->x_dim2; j++)
						if(randy <= hmm->place_M_cum[BIT(data->x_vec[t-1] ,0)][j][t-1])
						{
							data->x_vec[t] = j;
							break;
						}
					// determine the copy number (discrete data, between 0 and 7, stored in bits 1-3) for alpha
					randy = randr();
					for(j = 0; j < hmm->x_dim; j++)
						if(randy <= hmm->M_cum[BITS(data->x_vec[t-1],1,3)][j])
						{
							data->x_vec[t] ^= (j << 1);
							break;
						}					
					// determine the genotype (binary data, stored in the lsb) for beta
					randy = randr();
					for(j = 0; j < hmm->x_dim2; j++)
						if(randy <= hmm->place_M_cum[BIT(data->x_vec[t-1],4)][j][t-1])
						{
							data->x_vec[t] ^= (j << 4);
							break;
						}
					// determine the copy number (discrete data, between 0 and 7, stored in bits 5-7) for beta
					randy = randr();
					for(j = 0; j < hmm->x_dim; j++)
						if(randy <= hmm->M_cum[BITS(data->x_vec[t-1],5,7)][j])
						{
							data->x_vec[t] ^= (j << 5);
							break;
						}
//					data->x_vec[t] = 999;
				}
				// Artificial 'uniparental disomy' make one chromosome zero and the other chromosome 2
				// for(t=200; t < data->seq_len; t++)
				//	data->x_vec[t] = (data->x_vec[t]&0xFFFFFF11)|0x4;
			}	
			else // here do 'regular' HMM
				for(t = 1; t < data->seq_len; t++)
				{
					randy = randr();
					for(j = 0; j < hmm->x_dim; j++)
						if(randy <= hmm->place_M_cum[data->x_vec[t-1]][j][t-1])
						{
							data->x_vec[t] = j;
							break;
						}
				}
		}
		else  // Do the usual, with just one transition matrix
		{
			for(t = 1; t < data->seq_len; t++)
			{
				randy = randr();
				for(j = 0; j < hmm->x_dim; j++)
					if(randy <= hmm->M_cum[data->x_vec[t-1]][j])
					{
						data->x_vec[t] = j;
						break;
					}
			}
		}
	}		
	else   // here missing data
	{
		for(t = 1; t < data->seq_len; t++)
		{
			// first decide when to put it ... 
			interval = Geometric(data->dont_miss_prob); // distribution starts from zero ! 					
			data->loc_vec[t] = data->loc_vec[t-1] + interval + 1; //set the loc
			data->loc_diff_vec[t] = interval;  // Note : No gap is denoted by diff=0
		
			randy = randr();
			for(j = 0; j < hmm->x_dim; j++)
				if(randy < hmm->M_POWERS_cum[interval][data->x_vec[t-1]][j])
				{
					data->x_vec[t] = j;
					break;
				}
		}
	}

//return 2; // forget about the y's

	// Now simulate the y's
	if(hmm->y_type == DISCRETE)	
		for(t = 0; t < data->seq_len; t++)
		{
			randy = randr();
			for(j = 0; j < hmm->y_dim; j++)
				if(randy < hmm->N_cum[data->x_vec[t]][j])
				{
					data->y_vec[t] = j;
					data->y_vec_int[t] = j;
					break;
				}
		}
	else  // Here Mixture of Gaussians
	{
		if(hmm->special_models_flag == SNPS_ONE_SAMPLE)
			for(t = 0; t < data->seq_len; t++)
			{
				// First randomize which of the mixtures to take
				randy = randr();
				A_copy = (1^BIT(data->x_vec[t],0)) * BITS(data->x_vec[t],1,3) + 
							(1^BIT(data->x_vec[t],4)) * BITS(data->x_vec[t],5,7);
				B_copy = BIT(data->x_vec[t],0) * BITS(data->x_vec[t],1,3) + 
							 BIT(data->x_vec[t],4) * BITS(data->x_vec[t],5,7);
				mix_choose = 0;
				for(j = 0; j < hmm->y_dim; j++)
					if(randy <= hmm->N_cum[B_copy][j])
					{
						mix_choose = j;
						data->mix_vec[t] = j; // update what was chosen !!! 
						break;
					}
				// Now randomize according to a gaussian distribution 
				data->y_vec[t] = Gaussian(hmm->MU[A_copy][mix_choose], hmm->SIGMA[A_copy][mix_choose]);
				data->y_vecB[t] = Gaussian(hmm->MU[B_copy][mix_choose], hmm->SIGMA[B_copy][mix_choose]);
			}
		else  // Do the usual, with just one transition matrix
			for(t = 0; t < data->seq_len; t++)
			{
				// First randomize which of the mixtures to take
				randy = randr();
				for(j = 0; j < hmm->y_dim; j++)
					if(randy <= hmm->N_cum[data->x_vec[t]][j])
					{
						mix_choose = j;
						data->mix_vec[t] = j; // update what was chosen !!! 
						break;
					}

				// Now randomize according to a gaussian distribution 
				data->y_vec[t] = Gaussian(hmm->MU[data->x_vec[t]][mix_choose], 
					hmm->SIGMA[data->x_vec[t]][mix_choose]);
			}
	}
	return 0;

}




// Compute an approximate Kullback-Leibler Divergance rate between the two HMMs 
// Note : We do not assume that the data has enough allocation, so we allocate inside !!!  
// Any data already stored is ruined
long ComputeKLDistance(hmm_model *hmm1, hmm_model *hmm2, hmm_data *data, long iters, double *KL_dist)
{
	long i, j, k;

	double y_vec_log_prob_out_hmm1, y_vec_log_prob_out_hmm2;  // prob. according to both models

	data->miss_data = 0; // No missing data here !!!! 


	// Now allocate to the data struct 
	data->y_vec = new double[data->seq_len];
	data->y_vecB = new double[data->seq_len];

	data->y_vec_int = new long[data->seq_len];
	data->x_vec = new long[data->seq_len];
	data->mix_vec = new long[data->seq_len];
	data->loc_vec = new long[data->seq_len];
	data->loc_diff_vec = new long[data->seq_len];

	// allocate neccessary variables
	double *y_cond_tabs[MAX_X_VALS];
	double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS];
	double *alpha[MAX_X_VALS];
	double *scale = new double[data->seq_len];   // scaling to avoid underflows
	double *scale_cum = new double[data->seq_len];   // cumulative sum of the scale
	double *scale_exp = new double[data->seq_len];   // exp of the scale (ratios)

	// First allocate memory 
	for(i = 0; i < (MAX(hmm1->x_dim, hmm2->x_dim)); i++)
	{
		y_cond_tabs[i] = new double[data->seq_len];
		alpha[i] = new double[data->seq_len];	
		for(j = 0; j < (MAX(hmm1->y_dim, hmm2->y_dim)); j++)
			y_mu_square_exp_vecs[i][j] = new double[data->seq_len];
	}


	*KL_dist = 0; // init to zero
	printf("X DIMS : HMM1 %ld HMM2 %ld\n", hmm1->x_dim, hmm2->x_dim); 
	

	for(k = 0; k < iters; k++)
	{


		// simulate a sequence according to 1st model
		SimulateSequenceFromModel(hmm1, data);

/*******/

		// Compute prob. according to first model 
		Compute_y_cond_tabs(hmm1,  data, y_cond_tabs, y_mu_square_exp_vecs); 
		forward( hmm1, data, y_cond_tabs, 
			  alpha, scale, scale_cum, scale_exp, &y_vec_log_prob_out_hmm1);
		printf("Iter : %ld LogProb1 : %lf  ", k, y_vec_log_prob_out_hmm1);

		// Compute prob. according to second model 
		Compute_y_cond_tabs(hmm2,  data, y_cond_tabs, y_mu_square_exp_vecs); 
		forward( hmm2, data, y_cond_tabs, 
			  alpha, scale, scale_cum, scale_exp, &y_vec_log_prob_out_hmm2);
/*******/	
		printf("LogProb2 : %lf  ", y_vec_log_prob_out_hmm2);

		*KL_dist += (y_vec_log_prob_out_hmm1-y_vec_log_prob_out_hmm2) / data->seq_len;

	}

	(*KL_dist) /= iters;  // divide to get the average

	// finally delete all allocated buffers
	delete data->y_vec;
	delete data->y_vecB;

	delete data->y_vec_int;  
	delete data->x_vec;
	delete data->mix_vec;
	delete data->loc_vec;
	delete data->loc_diff_vec;
	delete scale;
	delete scale_cum;
	delete scale_exp;
	for(i = 0; i < (MAX(hmm1->x_dim, hmm2->x_dim)); i++)
	{
		delete y_cond_tabs[i];
		delete alpha[i];
		for(j = 0; j < (MAX(hmm1->y_dim, hmm2->y_dim)); j++)
			delete y_mu_square_exp_vecs[i][j];
	}


	return 0;

}





// Compute an approximate Kullback-Leibler Divergance rate between the two HMMs SNPs
// Note : We do not assume that the data has enough allocation, so we allocate inside !!!  
// Any data already stored is ruined
long ComputeKLDistanceSNPs(hmm_model *hmm1, hmm_model *hmm2, hmm_data *data, long iters, double *KL_dist)
{
	long j, k, i1, i2, i_geno, total_i_index;

	double y_vec_log_prob_out_hmm1, y_vec_log_prob_out_hmm2;  // prob. according to both models

	data->miss_data = 0; // No missing data here !!!! 


	// Now allocate to the data struct 
	data->y_vec = new double[ data->seq_len];
	data->y_vecB = new double[data->seq_len];
	data->x_vec = new long[data->seq_len];
	data->mix_vec = new long[data->seq_len];
	data->loc_vec = new long[data->seq_len];
	data->loc_diff_vec = new long[data->seq_len];

	// allocate neccessary variables
	double *y_cond_tabs[MAX_2D_X_VALS];
	double *y_cond_tabsB[MAX_2D_X_VALS];
	double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS];
	double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS];
	double *alpha[36];
	double *scale = new double[data->seq_len];   // scaling to avoid underflows
	double *scale_cum = new double[data->seq_len];   // cumulative sum of the scale
	double *scale_exp = new double[data->seq_len];   // exp of the scale (ratios)


	// First allocate memory 

	for(i_geno = 0; i_geno < (hmm1->x_dim2*hmm1->x_dim2); i_geno++)		// state at time t-1
		for(i1 = 0; i1 < hmm1->x_dim; i1++)		// state at time t-1
			for(i2 = 0; i2 < hmm1->x_dim; i2++)		// state at time t-1
			{
				total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 
				alpha[total_i_index] = new double[data->seq_len];	
				if(total_i_index <= 4)
				{
					y_cond_tabs[total_i_index] = new double[data->seq_len];
					y_cond_tabsB[total_i_index] = new double[data->seq_len];
					for(j = 0; j < (MAX(hmm1->y_dim, hmm2->y_dim)); j++)
					{
						y_mu_square_exp_vecs[total_i_index][j] = new double[data->seq_len];
						y_mu_square_exp_vecsB[total_i_index][j] = new double[data->seq_len];
					}
				}
			}


	*KL_dist = 0; // init to zero
	printf("X DIMS : HMM1 %ld HMM2 %ld\n", hmm1->x_dim, hmm2->x_dim); 
	
//////////////////////////////////////////////////////////////////////////////  Phase 0 ////////////////
/****************************************************************************/


	for(k = 0; k < iters; k++)
	{
		printf("DO ITER %ld\n", k); 
		// simulate a sequence according to 1st model
		SimulateSequenceFromModel(hmm1, data);
		// Compute prob. according to first model 
		Compute_y_cond_tabsSNPs(hmm1,  data, y_cond_tabs, y_cond_tabsB,
				 y_mu_square_exp_vecs, y_mu_square_exp_vecsB); // last repeat (dummy)		
		forwardSNPs( hmm1, data, y_cond_tabs, y_cond_tabsB, 
			  alpha, scale, scale_cum, scale_exp, &y_vec_log_prob_out_hmm1);
		// Compute prob. according to second model 
		Compute_y_cond_tabsSNPs(hmm2,  data, y_cond_tabs, y_cond_tabsB, 
			y_mu_square_exp_vecs, y_mu_square_exp_vecsB); // last repeat (dummy)		
		forwardSNPs( hmm2, data, y_cond_tabs, y_cond_tabsB,
			  alpha, scale, scale_cum, scale_exp, &y_vec_log_prob_out_hmm2);
		*KL_dist += (y_vec_log_prob_out_hmm1-y_vec_log_prob_out_hmm2) / data->seq_len;
	}

 /****************************************************************************/
//////////////////////////////////////////////////////////////////////////////  End Phase 0 ////////////////

	(*KL_dist) /= iters;  // divide to get the average

	// finally delete all allocated buffers
	delete data->y_vec;
	delete data->y_vecB;  
	delete data->x_vec;
	delete data->mix_vec;
	delete data->loc_vec;
	delete data->loc_diff_vec;
	delete scale;
	delete scale_cum;
	delete scale_exp;

	for(i_geno = 0; i_geno < (hmm1->x_dim2*hmm1->x_dim2); i_geno++)		// state at time t-1
		for(i1 = 0; i1 < hmm1->x_dim; i1++)		// state at time t-1
			for(i2 = 0; i2 < hmm1->x_dim; i2++)		// state at time t-1
			{
				total_i_index = multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]; 
				delete alpha[total_i_index];
				if(total_i_index <= 4)
				{
					delete y_cond_tabs[total_i_index];
					delete y_cond_tabsB[total_i_index];
					for(j = 0; j < (MAX(hmm1->y_dim, hmm2->y_dim)); j++)
						delete y_mu_square_exp_vecs[total_i_index][j];
				}
			}


	return 0;

}


// This is only done for the continunuous case : Here we replace the tabs of N by computing
// Mixture of Gaussians coefficients !! 
// Here b_cond_tabs[i][t] = Pr(Y = y_vec[t] | X = i = \sum_{j}  N[i][j] * Gaussian_{i,j} (y_vec[t]) )
// Note : If we want tommorow a different continuous distribution, we chagne only (?) here !!! 
// We compute here also for discrete distributions to gain (?) performance !!! 
long Compute_b_cond_tabs(hmm_model *hmm, double *y_vec, long *y_vec_int, long seq_len, 
					double *b_cond_tabs[MAX_X_VALS], double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS])
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
						y_mu_square_exp_vecs[i][j][t] =  ONE_OVER_2SQRT_PI * hmm->N[i][j] * hmm->SIGMA_inv[i][j] *
						exp( -(y_vec[t] - hmm->MU[i][j])*(y_vec[t] - hmm->MU[i][j]) * 0.5 *hmm->SIGMA_inv[i][j]*hmm->SIGMA_inv[i][j] ); 
		else   // here different gaussian for every place 
			for(i = 0; i < hmm->x_dim; i++)
				for(j = 0; j < hmm->y_dim; j++)
					for(t = 0; t < seq_len; t++)
						y_mu_square_exp_vecs[i][j][t] =  ONE_OVER_2SQRT_PI * hmm->N[i][j] * 
						exp( -(y_vec[t] - hmm->place_gaussian_mu[i][j][t])*(y_vec[t] - hmm->place_gaussian_mu[i][j][t]) / 
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
					exp( -(y_vec[t] - hmm->MU[i][j])*(y_vec[t] - hmm->MU[i][j]) / (2.0*hmm->SIGMA[i][j]*hmm->SIGMA[i][j]) ) / 
					hmm->SIGMA[i][j];

*/
						b_cond_tabs[i][t] += 
					//	hmm->N[i][j] * 
						y_mu_square_exp_vecs[i][j][t];/// *
					//	exp( -y_mu_square_exp_vecs[i][j][t] * 
					//	hmm->SIGMA_inv[i][j]*hmm->SIGMA_inv[i][j]*0.5 ) * 
					//	hmm->SIGMA_inv[i][j];
					

					}						
	
				//////	b_cond_tabs[i][t] *= ONE_OVER_2SQRT_PI;
		//					exp( hmm->N_log[i][j] + ONE_OVER_2SQRT_PI_LOG + 
		//					 (y_vec[t] - hmm->MU[i][j])*(y_vec[t] - hmm->MU[i][j]) / (2.0*hmm->SIGMA[i][j]*hmm->SIGMA[i][j])  -
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
