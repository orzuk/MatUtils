#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "general.h"
#include "hmm_chrom_funcs.h"
#include "hapmap_funcs.h"



// Test several of the functions
main()
{
	
	long i, j, t;
	
	long x_dim = 3; 
	long y_dim = 1; // We want only one gaussian !!!
	long y_type = CONTINUOUS; // ; // CONTINUOUS; // DISCRETE; // CONTINUOUS;   // Gaussian
	
	long place_flag = 0; // 1 different gaussian at every place 
	long special_models_flag = 0; // 1 means the special SNPs model
	long use_bounds = 0; // do not constrain the parameters
	long miss_data = 0;  // is there data missing ??? 
	long fold_change_flag = 1; // No fold change for now !!! 
	double dont_miss_prob = 0.95; // prob. of getting the word
	long seq_len = 25700; 

	long copy_data = 0;  // flag saying if to copy the data

	double *y_vec = new double[seq_len];
	long *y_vec_int = new long[seq_len];
	long *x_vec = new long[seq_len];
	long *mix_vec = new long[seq_len]; // which mixture to choose
	long *trained_x_vec = new long[seq_len];
	long *loc_vec = new long[seq_len];   // location of the Y's and X's in case of missing data ..

	double *y_cond_tabs[MAX_X_VALS];
	double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS];

	

	// Here we use a model where at every point we have a differetn gaussian

	hmm_model  my_hmm;
	hmm_model  trained_hmm;
	
	hmm_data   my_data; // a structure containing all the data


// Compute prob. of the output sequence
	double y_vec_log_prob;

	double *alpha[MAX_X_VALS];
	double *beta[MAX_X_VALS];
	double *gamma[MAX_X_VALS];
	double *gammagamma[MAX_X_VALS][MAX_Y_VALS];
	double *scale = new double[seq_len];   // scaling to avoid underflows
	double *scale_cum = new double[seq_len];   // cumulative sum of the scale
	double *scale_exp = new double[seq_len];   // exp of the scale (ratios)
	long *count_occur[MAX_X_VALS][MAX_Y_VALS];


//	long timtim;   
	srand(12200); //   
//	srand(long(time(&timtim)));  // randomize seed 



	// Try reading HAPMAP data
//	long num_lines;
	FILE *genotype_f;
//	char *genotype_p = "E:\Research\HMM_Chromosome\SNP\HAPMAP\Genotype_data\WS_FTP.txt"; // genotypes_chr22_CEU_r21_nr_fwd.txt";	
	char *genotype_p = "genotypes_chr22_CEU_r21_nr_fwd.txt";	


	snp_hapmap_data snp_hapmap;



	long i_geno, ii_geno, i1, i2, total_i_index, i_geno1, i_geno2;
	long ctr=0;
	for(i_geno = 0; i_geno < 4; i_geno++)		// state at time t-1
		for(i1 = 0; i1 < 3; i1++)		// state at time t-1
			for(i2 = 0; i2 < 3; i2++)		// state at time t-1	
			{
				// Seperate according to whether the source or destination are zero ..
				total_i_index = ctr++; //// BIT(i_geno,0) + (BIT(i_geno,1) << 4) + (i1 << 1) + (i2 << 5);
//				total_j_index = BIT(j_geno,0) + (BIT(j_geno,1) << 4) + (j1 << 1) + (j2 << 5);
//				
//				multi_dim_total_index_tab[BIT(j_geno,0)][BIT(j_geno,1)][j1][j2]
				
					printf("%ld, |  %ld \n", total_i_index, multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)] );
			}
	printf("\n\n");	
	for(i_geno = 0; i_geno < 2; i_geno++)		// state at time t-1
		for(ii_geno = 0; ii_geno < 2; ii_geno++)		// state at time t-1
			for(i1 = 0; i1 < 3; i1++)		// state at time t-1
				for(i2 = 0; i2 < 3; i2++)		// state at time t-1	
			{
				// Seperate according to whether the source or destination are zero ..
				total_i_index = BIT(i_geno,0) + (BIT(ii_geno,0) << 4) + (i1 << 1) + (i2 << 5);
//				total_j_index = BIT(j_geno,0) + (BIT(j_geno,1) << 4) + (j1 << 1) + (j2 << 5);
//				
//				multi_dim_total_index_tab[j1][BIT(j_geno,0)][j2][BIT(j_geno,1)]
				
					printf("%ld, ", total_i_index);
			}

	long A_copy, B_copy;
	printf("\nNow the A COPY\n");
	for(i_geno = 0; i_geno < 4; i_geno++)		// state at time t-1
	{
		for(i1 = 0; i1 < 3; i1++)		// state at time t-1
			for(i2 = 0; i2 < 3; i2++)		// state at time t-1	
			{
				// Seperate according to whether the source or destination are zero ..
					B_copy = (1^BIT(i_geno,0)) * i1 + (1^BIT(i_geno,1)) * i2; // B copy number
					A_copy = BIT(i_geno,0) * i1 + BIT(i_geno,1) * i2;	// A copy numnber
					printf("%ld ", A_copy );
			}
			printf("\n");
	}
	printf("\n A COPY 1 DIM:\n");


	for(total_i_index = 0; total_i_index < 36; total_i_index++)
	{
		for(i_geno = 0; i_geno < 4; i_geno++)		// state at time t-1
			for(i1 = 0; i1 < 3; i1++)		// state at time t-1
				for(i2 = 0; i2 < 3; i2++)		// state at time t-1	
					if( multi_dim_total_index_tab[i1][BIT(i_geno,0)][i2][BIT(i_geno,1)]  == total_i_index )
					{

//		i_geno = BIT(total_i_index,0) + 2*BIT(total_i_index,4);
//		i1 = BITS(total_i_index, 1, 3);
//		i2 = BITS(total_i_index, 5, 3);

						A_copy = (1^BIT(i_geno,0)) * i1 + (1^BIT(i_geno,1)) * i2; // B copy number
						B_copy = BIT(i_geno,0) * i1 + BIT(i_geno,1) * i2;	// A copy numnber
						printf("%ld, ", A_copy);

					}
	}



	printf("\n Second Try A COPY 1 DIM:\n");
	for(i1 = 0; i1 < 3; i1++)		// state at time t-1
		for(i_geno1 = 0; i_geno1 < 2; i_geno1++)		// state at time t-1
			for(i2 = 0; i2 < 3; i2++)		// state at time t-1	
				for(i_geno2 = 0; i_geno2 < 2; i_geno2++)		// state at time t-1
				{

					A_copy = (1^i_geno1) * i1 + (1^i_geno2) * i2; // B copy number
					B_copy = i_geno1 * i1 + i_geno2 * i2;	// A copy numnber
					printf("%ld, ", A_copy);

				}

	printf("\n Second Try B COPY 1 DIM:\n");
	for(i1 = 0; i1 < 3; i1++)		// state at time t-1
		for(i_geno1 = 0; i_geno1 < 2; i_geno1++)		// state at time t-1
			for(i2 = 0; i2 < 3; i2++)		// state at time t-1	
				for(i_geno2 = 0; i_geno2 < 2; i_geno2++)		// state at time t-1
				{

					A_copy = (1^i_geno1) * i1 + (1^i_geno2) * i2; // B copy number
					B_copy = i_geno1 * i1 + i_geno2 * i2;	// A copy numnber
					printf("%ld, ", B_copy);

				}



	exit(99);

	snp_hapmap.num_lines = CountFileLines(genotype_f, genotype_p);
	printf("File contains %d lines\n", snp_hapmap.num_lines);


	long A=0x12345678;
	printf("A pop = %ld\n", PopCount(A));


	// allocate memory
	snp_hapmap.names = new char*[snp_hapmap.num_lines]; 
	snp_hapmap.data = new long*[snp_hapmap.num_lines]; 
	snp_hapmap.bases = new long[snp_hapmap.num_lines];
	snp_hapmap.freqs = new double[snp_hapmap.num_lines];
	snp_hapmap.hetros = new double[snp_hapmap.num_lines];
	snp_hapmap.bad_calls = new long*[snp_hapmap.num_lines];
	snp_hapmap.chrom_locs = new long[snp_hapmap.num_lines];


	for(i=0;i<snp_hapmap.num_lines;i++)
	{
		snp_hapmap.names[i] = new char[20];  // way more than maximal SNP name length
		snp_hapmap.data[i] = new long[6]; // 6 words gives 192 which is more than 180 bits needed (the two most are not used)
		snp_hapmap.bad_calls[i] = new long[6]; // 6 words gives 192 which is more than 180 bits needed (the two most are not used)
	}
		
	ReadGenotypeFile(genotype_f, genotype_p, &snp_hapmap); 


	// Printf the first part of the data, see if its the same
	for(i=0;i<50;i++)
	{
		for(j=0;j<6;j++)
			printf("%lf ", double(snp_hapmap.data[i][j]));
		printf("\n");
	}
	printf("num lines outside is %d\n", snp_hapmap.num_lines);

	// free memory
	for(i=0;i<snp_hapmap.num_lines;i++)
	{
		delete snp_hapmap.names[i];
		delete snp_hapmap.data[i];
		delete snp_hapmap.bad_calls[i];
	}		
	delete snp_hapmap.names; 
	delete snp_hapmap.data; 
	delete snp_hapmap.bases;
	delete snp_hapmap.freqs;
	delete snp_hapmap.hetros;
	delete snp_hapmap.bad_calls;
	delete snp_hapmap.chrom_locs;

	exit(99);








	printf("INITILIZING PROGRAM !\n");
	// Init the model
	my_hmm.cum_flag = 1; my_hmm.log_flag = 1;
	trained_hmm.cum_flag = 1; trained_hmm.log_flag = 1;


	InitilizeData(&my_data, seq_len, y_type, miss_data, dont_miss_prob, copy_data, 
					y_vec, y_vec_int, loc_vec);
	InitilizeModel(&my_hmm, x_dim, y_dim, y_type, use_bounds, place_flag, special_models_flag, miss_data, fold_change_flag);

	// a 'trick' to decieve the permuting until the arrays of mu are ready !!! 
	long remember_place = 0; 
	if(my_hmm.place_flag)	
	{
		my_hmm.place_flag = 0; 
		remember_place = 1; 
	}
	printf("Before permuting...\n");

	PermuteModelIncreasingMean(&my_hmm, &my_data); // assume we have a 'sorted' model 
	printf("After permuting...\n");
	
	my_hmm.place_flag = remember_place;  // bring it back !!! 

	PrintModel(&my_hmm);

	printf("Allocating ...\n");

	// First allocate memory 
	for(i = 0; i < my_hmm.x_dim; i++)
	{
	//	y_cond_tabs[i] = new double[seq_len];
		alpha[i] = new double[my_data.seq_len];
		beta[i] = new double[my_data.seq_len];
		gamma[i] = new double[my_data.seq_len];

		for(j = 0; j < my_hmm.y_dim; j++)
		{
			gammagamma[i][j] = new double[my_data.seq_len];
			count_occur[i][j] = new long[my_data.seq_len];
		}
	
	}

	printf("data allocating seq_len is %ld my_data_seq_len is %ld\n", seq_len, my_data.seq_len);

	// Now allocate to the data struct 
	my_data.y_vec = new double[ my_data.seq_len];
	my_data.y_vec_int = new long[my_data.seq_len];
	my_data.x_vec = new long[my_data.seq_len];
	my_data.mix_vec = new long[my_data.seq_len];
	my_data.loc_vec = new long[my_data.seq_len];
	my_data.loc_diff_vec = new long[my_data.seq_len];

	// Now the data mu and sigma for different gaussians
	if(my_hmm.place_flag) 
	{
		printf("alolo cc to zero\n");
		// allocate memory for mu and sigma 
		for(i = 0; i < my_hmm.x_dim; i++)
			for(j = 0; j < my_hmm.y_dim; j++)
			{
				my_hmm.place_gaussian_mu[i][j] = new double[my_data.seq_len];
				my_hmm.place_gaussian_sigma[i][j] = new double[my_data.seq_len];
			}


				printf("init to zero\n");
		// init to zero
		for(i = 0; i < my_hmm.x_dim; i++)
			for(j = 0; j < my_hmm.y_dim; j++)
				for(t = 0; t < my_data.seq_len; t++)
				{
					my_hmm.place_gaussian_mu[i][j][t] = 0;
					my_hmm.place_gaussian_sigma[i][j][t] = 0;
					count_occur[i][j][t] = 0; 
				}
		long sims;



		// Now try to determine the values 
		printf("NOW SIMULATING :\n");
		for(sims = 0; sims < 500; sims++)  // simulate 20 times for determining the means ..
		{
			SimulateSequenceFromModel(&my_hmm, &my_data);

			// add the value to mu 
			for(t = 0; t < my_data.seq_len; t++)
				my_hmm.place_gaussian_mu[my_data.x_vec[t]][mix_vec[t]][t] += my_data.y_vec[t];

			// add the value to sigma
			for(t = 0; t < my_data.seq_len; t++)
				my_hmm.place_gaussian_sigma[my_data.x_vec[t]][mix_vec[t]][t] += my_data.y_vec[t]*my_data.y_vec[t];				
	
			// count the number of occourences 
			for(t = 0; t < my_data.seq_len; t++)
				count_occur[my_data.x_vec[t]][mix_vec[t]][t] ++;				
		}

		// Now divide by number of occurences
		for(i = 0; i < my_hmm.x_dim; i++)
			for(j = 0; j < my_hmm.y_dim; j++)
			{
				for(t = 0; t < my_data.seq_len; t++)
				{
					if(count_occur[i][j][t] > 0)
					{
						my_hmm.place_gaussian_mu[i][j][t] /= double(count_occur[i][j][t]);
						my_hmm.place_gaussian_sigma[i][j][t] /= double(count_occur[i][j][t]);
					}
					else // we haven't encountered this count ... 
					{
						my_hmm.place_gaussian_mu[i][j][t] = 0.0; 
						my_hmm.place_gaussian_sigma[i][j][t] = 1.0; 
					}

				}
			}
	
		// Now subtract : sig = sqrt( E X^2 - (E X)^2 )
		for(i = 0; i < my_hmm.x_dim; i++)
			for(j = 0; j < my_hmm.y_dim; j++)
				for(t = 0; t < my_data.seq_len; t++)
					my_hmm.place_gaussian_sigma[i][j][t] = MAX(
					sqrt(my_hmm.place_gaussian_sigma[i][j][t] - 
					my_hmm.place_gaussian_mu[i][j][t] * my_hmm.place_gaussian_mu[i][j][t]), EPSILON); 


		// print some mus and sigmas
		for(i = 0; i < my_hmm.x_dim; i++)
			for(j = 0; j < my_hmm.y_dim; j++)
			{
				printf(" (I J) = (%ld %ld)\n", i, j);  
				for(t = 0; t < 20; t++)
					printf("%lf  |  %lf\n", my_hmm.place_gaussian_mu[i][j][t], 
					my_hmm.place_gaussian_sigma[i][j][t]);
			}



		// Now do again this increasing order stuff !!!!! 
		PermuteModelIncreasingMean(&my_hmm, &my_data); // assume we have a 'sorted' model 

	}



	printf("Now Simulating\n");
	// simulate a sequence   
	SimulateSequenceFromModel(&my_hmm, &my_data);
	printf("finished sim\n");
	// Copy to integers 
//	for(t = 0; t < my_data.seq_len; t++)
//		y_vec_int[t] = long(y_vec[t]);

	// allocate memory
	for(i = 0; i < my_hmm.x_dim; i++)
		y_cond_tabs[i] = new double[my_data.seq_len];
	for(i = 0; i < my_hmm.x_dim; i++)
		for(j = 0; j < my_hmm.y_dim; j++)
			y_mu_square_exp_vecs[i][j] = new double[my_data.seq_len];





	// first need to compute the helping b_tabs
	Compute_y_cond_tabs(&my_hmm, &my_data, y_cond_tabs, /*y_cond_tabs, y_mu_square_exp_vecs,*/ y_mu_square_exp_vecs);
	printf("finished b_cond\n");

	double x_loglike, x_place_loglike, y_loglike;
	// Determine best X vec generating the Y's vec ..
	Viterbi(&my_hmm, &my_data, y_cond_tabs, y_cond_tabs, trained_x_vec, &x_loglike, &x_place_loglike, &y_loglike); // dummy repeat




//#define PRINT2
#ifdef PRINT2
	// Print results ...
	printf("The original X : \n");	print_vec(x_vec, my_data.seq_len);
//	if(my_hmm.y_type == DISCRETE)
//	{
///		for(t = 0; t < my_data.seq_len; t++)
//			y_vec_int[t] = long(y_vec[t]);
//		printf("The observed Y : \n");	print_vec(y_vec_int, my_data.seq_len);
//	}
//	else
//		printf("The observed Y : \n");	print_double_vec(y_vec, my_data.seq_len);

	printf("The best-path X : \n");	print_vec(trained_x_vec, my_data.seq_len);
#endif

	// Count match
	long match_count = 0;
	for(i = 0; i < my_data.seq_len; i++)
		if(my_data.x_vec[i] == trained_x_vec[i])
			match_count++;

	printf("Found %ld Matches out of %ld : Gives %lf %% !!\n", match_count, my_data.seq_len, 100.0*double(match_count)/double(my_data.seq_len));





	

	// first need to compute the helping b_tabs
	Compute_y_cond_tabs(&my_hmm, &my_data, y_cond_tabs, /*y_cond_tabs, y_mu_square_exp_vecs,*/ y_mu_square_exp_vecs);


	// Print the b's
	/*
	printf("The B_COND_TABS : \n");
	for(t = 0; t < seq_len; t++)
	{
		for(i = 0; i < x_dim; i++)
			printf("%lf ", y_cond_tabs[i][t]);
		printf("\n");
	}
	printf("-----------\n\n");
	*/


	forward(&my_hmm, &my_data, y_cond_tabs,
			  alpha, scale, scale_cum, scale_exp, &y_vec_log_prob);
	printf("Forward Sequence log-prob is : %lf \n", y_vec_log_prob);
	
	backward(&my_hmm, &my_data, scale_exp, y_cond_tabs,
			  beta);
//	printf("Backward Sequence log-prob is : %lf \n", y_vec_log_prob);

	ComputeMarginalGamma(&my_hmm, &my_data, y_mu_square_exp_vecs,alpha, beta, scale_exp, y_cond_tabs, 
						 gamma, gammagamma);




	printf("Finished Gamma !!!\n");
#ifdef PRINT2
	printf("Conditional probs : \n");
	for(i = 0; i < my_hmm.x_dim; i++)
	{
		print_double_vec(y_cond_tabs[i], my_data.seq_len);
		printf("\n-----------------------\n");
	}

	printf("Marginal probs : \n");
	for(i = 0; i < my_hmm.x_dim; i++)
		print_double_vec(gamma[i], my_data.seq_len);
#endif

//#ifdef LATER	

	long max_iters = 80; // number of iterations 
	double tolerance = 0.00001;
	long num_starts = 8;  // different random places to start from
	double best_model_score;
	

	// Try to learn parameters back
	InitilizeModel(&trained_hmm, x_dim, y_dim, y_type, use_bounds, place_flag, special_models_flag, miss_data, fold_change_flag);


	if(trained_hmm.place_flag) 
	{
		// allocate memory for mu and sigma 
		for(i = 0; i < my_hmm.x_dim; i++)
			for(j = 0; j < my_hmm.y_dim; j++)
			{
				trained_hmm.place_gaussian_mu[i][j] = new double[my_data.seq_len];
				trained_hmm.place_gaussian_sigma[i][j] = new double[my_data.seq_len];
			}

		// Copy the different mu and sigma
		CopyModel(&my_hmm, &trained_hmm, &my_data); 

		// init again to randomize all other parameters
		InitilizeModel(&trained_hmm, x_dim, y_dim, y_type, use_bounds, place_flag, special_models_flag, miss_data, fold_change_flag);

	}
 

	printf("START EM Outside !!!\n");
	clock_t EM_start_time = clock();

	// location vec denotes the locations of the observations. We assume that between them we do not see anything !!! 
	TrainModelEM(&trained_hmm, &my_data, max_iters, num_starts, tolerance, &best_model_score);


	clock_t EM_end_time = clock();
	printf("END EM  !!! Score : %lf\n", best_model_score);

	printf("THE ORIGINAL MODEL : (log-score %lf ) \n", y_vec_log_prob);
	PrintModel(&my_hmm);

	Compute_y_cond_tabs(&trained_hmm, &my_data,  
		y_cond_tabs, /*y_cond_tabs, y_mu_square_exp_vecs,*/ y_mu_square_exp_vecs);
	forward(&trained_hmm, &my_data, y_cond_tabs, 
			  alpha, scale, scale_cum, scale_exp, &y_vec_log_prob);
	printf("THE FINAL TRAINED MODEL : (log-score %lf ) \n", y_vec_log_prob);
	PrintModel(&trained_hmm);

	printf("\n-------\nComparison between TRUE and TRAINED model :\n");
	// Now compare the learned model to the original one !
	CompareTwoModels(&my_hmm, &trained_hmm);

	
	// Now compare CORRECT model and a random model
	InitilizeModel(&trained_hmm, my_hmm.x_dim, my_hmm.y_dim, my_hmm.y_type, my_hmm.use_bounds, 
		my_hmm.place_flag, my_hmm.special_models_flag, my_hmm.miss_data, my_hmm.fold_change_flag);
	
	printf("\n---------\n Now comparing TRUE model to RANDOM model : \n-------\n");
	CompareTwoModels(&my_hmm, &trained_hmm);



	// Now compare two random models and see what we get 
	InitilizeModel(&my_hmm, my_hmm.x_dim, my_hmm.y_dim, my_hmm.y_type, my_hmm.use_bounds, 
		my_hmm.place_flag, my_hmm.special_models_flag, my_hmm.miss_data, my_hmm.fold_change_flag);
	InitilizeModel(&trained_hmm, my_hmm.x_dim, my_hmm.y_dim, my_hmm.y_type, my_hmm.use_bounds, 
		my_hmm.place_flag, my_hmm.special_models_flag, my_hmm.miss_data, my_hmm.fold_change_flag);
	
	printf("\n---------\n Now comparing two RANDOM models : \n-------\n");
	CompareTwoModels(&my_hmm, &trained_hmm);



	printf("Total EM Time (sec.) : %lf\n", double(EM_end_time-EM_start_time) / double(CLOCKS_PER_SEC));
//#endif


// #define CANCER
#ifdef CANCER

	printf("\n\n--------------\nNow Do Cancer ! ! ! \n-------------------\n\n");
	max_iters = 12; num_starts = 5;

/*	double cancer_exp_vec[316] = {1.234567 ,1.234567 ,8.625958 ,1.234567 ,4.447346 ,6.714049 ,1.234567 ,1.234567 ,
6.749111 ,7.144960 ,4.115780 ,1.234567 ,5.831882 ,7.051682 ,4.628887 ,6.337003 ,
6.234411 ,1.234567 ,5.777962 ,6.592496 ,1.234567 ,1.131402 ,4.207673 ,1.234567 ,
4.598146 ,5.247550 ,1.234567 ,8.231536 ,7.203926 ,1.234567 ,8.159804 ,5.473530 ,
3.990834 ,1.234567 ,6.797829 ,5.539694 ,1.234567 ,6.985827 ,1.234567 ,4.708629 ,
6.440787 ,1.234567 ,5.131672 ,5.027165 ,6.133181 ,7.011935 ,5.180097 ,8.387267 ,
1.234567 ,4.817859 ,6.233234 ,4.601162 ,4.583947 ,5.335613 ,1.234567 ,5.757955 ,
6.064250 ,1.234567 ,1.234567 ,3.968403 ,3.732896 ,1.234567 ,4.395683 ,6.567656 ,
1.234567 ,4.757891 ,5.819192 ,5.747799 ,4.833898 ,6.897099 ,1.234567 ,1.234567 ,
7.078932 ,6.351758 ,10.240231 ,7.775696 ,6.606785 ,6.675319 ,7.342909 ,1.234567 ,
5.584248 ,1.234567 ,8.519550 ,6.791558 ,6.932936 ,7.361312 ,6.886532 ,7.107425 ,
5.652489 ,1.234567 ,6.215008 ,1.234567 ,1.234567 ,1.234567 ,1.234567 ,5.842094 ,
1.234567 ,4.127134 ,6.873681 ,3.640214 ,1.234567 ,7.014545 ,7.206674 ,1.234567 ,
5.642262 ,1.234567 ,6.841829 ,8.135933 ,1.234567 ,7.640700 ,1.234567 ,7.175949 ,
6.389233 ,2.778819 ,7.234971 ,3.020425 ,4.448516 ,4.618086 ,6.380631 ,5.060060 ,
6.980820 ,4.879007 ,6.731734 ,6.133398 ,5.502482 ,1.234567 ,6.662749 ,1.234567 ,
1.234567 ,7.393263 ,5.585749 ,6.299501 ,6.860244 ,5.557600 ,1.234567 ,7.841454 ,
6.029483 ,6.235782 ,6.475279 ,6.527080 ,7.580138 ,8.990628 ,1.234567 ,6.451891 ,
6.128614 ,6.215807 ,7.344848 ,1.234567 ,6.177114 ,9.072124 ,10.974953 ,1.234567 ,
5.593596 ,1.234567 ,5.409411 ,6.949665 ,8.177825 ,3.958907 ,6.159095 ,5.545568 ,
6.647299 ,6.280396 ,6.779240 ,5.479805 ,5.844124 ,1.234567 ,1.234567 ,4.742320 ,
1.234567 ,8.928442 ,6.634108 ,7.021620 ,6.798387 ,4.496471 ,1.234567 ,3.813307 ,
5.129307 ,4.908233 ,3.020425 ,1.234567 ,4.427239 ,5.879135 ,1.234567 ,6.260346 ,
6.832816 ,6.757281 ,6.120737 ,5.695078 ,1.234567 ,4.055257 ,6.843217 ,1.234567 ,
6.970918 ,1.234567 ,1.234567 ,9.383360 ,7.133535 ,5.262172 ,5.374352 ,6.793690 ,
5.348059 ,1.234567 ,1.234567 ,8.850704 ,8.950494 ,6.805058 ,1.234567 ,6.756002 ,
7.032889 ,7.730921 ,6.193180 ,6.352455 ,8.636131 ,5.796969 ,6.413623 ,5.270946 ,
6.650279 ,7.579372 ,5.662613 ,4.666265 ,7.702917 ,6.285998 ,6.009550 ,7.639257 ,
1.234567 ,3.295837 ,7.024827 ,10.393545 ,6.323642 ,6.782985 ,1.234567 ,6.940706 ,
8.318303 ,1.234567 ,5.283711 ,5.963579 ,7.430292 ,6.002652 ,5.565286 ,6.683611 ,
4.802380 ,3.580737 ,7.473694 ,1.234567 ,7.072676 ,4.763028 ,4.997888 ,4.779123 ,
1.234567 ,7.240936 ,6.135781 ,7.815207 ,6.828171 ,7.816175 ,6.561455 ,4.650144 ,
6.546068 ,3.795489 ,1.234567 ,6.975507 ,6.856778 ,9.286458 ,3.962716 ,1.234567 ,
6.078559 ,5.335131 ,4.670958 ,5.197391 ,9.325516 ,6.900731 ,6.295635 ,5.982172 ,
5.329331 ,7.347686 ,7.266758 ,7.091825 ,4.891852 ,5.931184 ,4.897093 ,1.234567 ,
5.741399 ,1.234567 ,8.256867 ,7.268502 ,6.563573 ,4.691348 ,5.778581 ,4.958640 ,
6.426974 ,5.274537 ,6.098074 ,6.229300 ,7.219788 ,1.234567 ,4.497585 ,6.044294 ,
4.850467 ,5.633360 ,3.230804 ,6.156555 ,1.234567 ,8.266781 ,6.976535 ,5.375741 ,
1.234567 ,5.719328 ,7.253187 ,5.788736 ,5.106551 ,4.070735 ,5.798790 ,1.234567 ,
5.351384 ,5.391352 ,3.284664 ,2.949688 };
*/






	double cancer_exp_vec[1104] = { 5.245444 ,0.000000 ,6.116334 ,8.564668 ,7.220593 ,0.000000 ,5.731073 ,7.134652 ,7.414573 ,7.392401 ,7.940370 ,8.831580 ,6.246882 ,
4.165114 ,6.487075 ,5.395898 ,7.752249 ,8.210179 ,5.730749 ,7.099037 ,7.613029 ,4.708629 ,9.525589 ,5.660875 ,0.000000 ,0.000000 ,
0.000000 ,0.000000 ,6.830010 ,7.301687 ,9.322534 ,0.000000 ,9.623906 ,7.334590 ,7.976939 ,6.414606 ,0.000000 ,5.085124 ,3.879500 ,
6.588238 ,0.000000 ,6.199291 ,8.487476 ,5.635504 ,7.783182 ,4.809742 ,7.361121 ,8.436980 ,5.041488 ,8.356132 ,4.332048 ,9.782652 ,
0.000000 ,0.000000 ,6.973075 ,6.454569 ,0.000000 ,5.785977 ,9.036380 ,6.425355 ,0.000000 ,0.000000 ,7.167269 ,7.659124 ,
5.087596 ,5.254888 ,6.399593 ,6.930299 ,6.033566 ,8.181665 ,5.026509 ,0.000000 ,5.521461 ,8.650079 ,6.355934 ,6.909354 ,
5.691710 ,4.082609 ,6.402249 ,6.424545 ,6.143971 ,0.000000 ,7.835184 ,0.000000 ,7.177401 ,6.149536 ,7.544068 ,0.000000 ,
3.749504 ,8.766379 ,0.000000 ,6.484177 ,6.711010 ,2.104134 ,6.761457 ,5.107762 ,4.891852 ,8.708953 ,6.847899 ,3.456317 ,
6.328472 ,5.674354 ,6.922545 ,6.096050 ,5.508578 ,8.179704 ,6.878944 ,6.562162 ,6.709060 ,0.000000 ,6.579112 ,7.025272 ,
0.000000 ,8.471589 ,4.457830 ,8.447007 ,8.169648 ,7.187279 ,5.234845 ,7.113142 ,6.495568 ,7.442024 ,6.827413 ,7.883069 ,
8.028129 ,6.078788 ,3.384390 ,7.945733 ,8.673889 ,9.317337 ,6.270988 ,5.912151 ,7.626180 ,6.860349 ,6.978680 ,7.395353 ,
5.563754 ,9.075173 ,6.192158 ,0.000000 ,7.818108 ,5.790266 ,0.000000 ,0.000000 ,3.953165 ,7.287629 ,7.162708 ,7.603399 ,
0.000000 ,0.000000 ,5.140493 ,6.051148 ,4.018183 ,5.286751 ,7.168734 ,7.372684 ,0.000000 ,6.988967 ,3.919991 ,7.925736 ,
10.643447 ,0.262364 ,6.162472 ,6.512785 ,9.434587 ,7.359277 ,7.243370 ,7.795688 ,5.398163 ,8.284025 ,5.201806 ,6.602588 ,
7.281386 ,2.525729 ,5.896879 ,0.000000 ,6.189905 ,0.000000 ,6.003887 ,0.000000 ,6.647429 ,5.499215 ,4.825109 ,7.371615 ,
5.864483 ,3.691376 ,9.981467 ,6.444608 ,5.654592 ,0.000000 ,0.000000 ,6.083360 ,6.345110 ,7.280491 ,4.133565 ,0.000000 ,
7.992201 ,8.119547 ,5.751620 ,5.846728 ,5.855072 ,4.898586 ,5.155024 ,6.755653 ,6.502490 ,6.287673 ,4.879767 ,7.371426 ,
6.710158 ,0.336472 ,8.061834 ,0.000000 ,3.711130 ,6.933033 ,6.622736 ,7.389070 ,6.013471 ,0.000000 ,3.549617 ,3.756538 ,
8.317229 ,7.343685 ,4.968423 ,0.000000 ,7.734121 ,8.815325 ,7.242726 ,2.674149 ,0.000000 ,6.670513 ,6.947841 ,6.407045 ,
8.940616 ,8.250803 ,0.000000 ,3.339322 ,7.029619 ,7.592215 ,8.313386 ,8.050321 ,0.000000 ,6.210198 ,2.995732 ,7.453852 ,
8.017934 ,8.632235 ,7.083388 ,7.614756 ,0.000000 ,6.429881 ,8.414163 ,0.000000 ,7.014724 ,7.675175 ,5.905634 ,3.698830 ,
7.588425 ,6.641313 ,0.000000 ,7.091492 ,7.973638 ,8.768434 ,0.000000 ,7.604795 ,4.793308 ,4.439116 ,0.000000 ,0.000000 ,
6.628835 ,7.527202 ,3.394508 ,7.106688 ,6.555499 ,8.667628 ,4.606170 ,8.358971 ,0.000000 ,8.984293 ,6.632397 ,6.888471 ,
5.718671 ,0.000000 ,8.719677 ,5.510198 ,7.687355 ,5.633002 ,3.891820 ,5.304796 ,7.632013 ,5.742683 ,0.000000 ,0.000000 ,
0.000000 ,7.209118 ,5.357529 ,9.910915 ,7.623495 ,5.150397 ,6.432136 ,7.320659 ,7.313953 ,6.839155 ,6.844176 ,6.882335 ,
6.604350 ,6.452049 ,9.093795 ,6.655183 ,6.472037 ,7.178317 ,3.908015 ,10.832781 ,7.043684 ,7.733333 ,7.211188 ,5.697429 ,
4.517431 ,5.771753 ,5.795450 ,9.882443 ,8.446084 ,7.108571 ,0.000000 ,6.867037 ,0.000000 ,0.000000 ,5.617135 ,5.591733 ,
5.040841 ,7.344525 ,8.850862 ,5.082025 ,5.967172 ,6.330256 ,8.455509 ,7.169196 ,4.711330 ,4.044804 ,5.011302 ,2.351375 ,
6.122931 ,7.736482 ,0.000000 ,8.463012 ,8.589942 ,0.000000 ,6.352978 ,7.087991 ,6.416078 ,6.116995 ,10.291325 ,6.835722 ,
7.387647 ,0.000000 ,6.692208 ,5.805737 ,2.734368 ,5.579730 ,6.626453 ,0.000000 ,7.363914 ,4.043051 ,7.873560 ,8.805195 ,
0.000000 ,5.590987 ,5.493473 ,6.686235 ,8.453145 ,0.000000 ,0.000000 ,5.739793 ,3.572346 ,3.906005 ,5.610570 ,5.554896 ,
4.286341 ,3.095578 ,6.072353 ,2.660260 ,6.093570 ,6.357669 ,6.765961 ,4.070735 ,8.851821 ,7.794823 ,7.108490 ,6.324000 ,
3.182212 ,7.915969 ,4.745801 ,0.000000 ,4.958640 ,7.662562 ,8.176195 ,6.460217 ,6.395595 ,7.348845 ,8.477287 ,3.725693 ,
0.000000 ,8.388268 ,0.000000 ,6.937799 ,4.728272 ,0.000000 ,0.000000 ,5.929855 ,0.000000 ,0.000000 ,5.892749 ,7.010943 ,
7.608672 ,7.305121 ,3.781914 ,7.271912 ,7.497318 ,5.353752 ,7.058414 ,5.610570 ,0.000000 ,4.206184 ,5.063860 ,5.558371 ,
0.000000 ,0.000000 ,6.313729 ,6.608270 ,1.609438 ,3.269569 ,0.000000 ,4.897093 ,2.631889 ,0.000000 ,7.980810 ,4.529368 ,
7.424821 ,8.176476 ,7.507086 ,5.035003 ,7.024649 ,6.776165 ,6.254982 ,5.085124 ,0.000000 ,5.439817 ,6.122273 ,0.000000 ,
5.218733 ,0.000000 ,6.958448 ,6.087456 ,3.100092 ,6.244361 ,10.369132 ,5.421420 ,6.354370 ,6.978866 ,5.718343 ,6.827738 ,
6.062389 ,5.864767 ,5.501258 ,5.070161 ,6.900932 ,4.922168 ,7.314619 ,5.709765 ,7.370797 ,6.226339 ,6.577722 ,10.315481 ,
10.646486 ,5.801212 ,5.879135 ,4.034241 ,5.981919 ,7.400071 ,6.840333 ,5.026509 ,6.421785 ,4.733563 ,0.000000 ,7.064759 ,
4.571613 ,4.822698 ,5.677438 ,6.489053 ,2.890372 ,2.091864 ,4.471639 ,6.782645 ,0.000000 ,5.199601 ,6.085410 ,10.406473 ,
8.502019 ,7.594432 ,3.443618 ,7.233311 ,5.301811 ,4.462454 ,7.530587 ,8.540597 ,8.114205 ,5.283711 ,6.006353 ,6.411982 ,
6.105909 ,6.092440 ,8.144708 ,9.015128 ,8.059150 ,7.758419 ,8.799405 ,0.000000 ,5.236442 ,2.912351 ,7.413729 ,8.541574 ,
0.000000 ,5.277094 ,0.000000 ,5.897154 ,6.512042 ,5.378514 ,6.355587 ,5.799396 ,7.418061 ,6.510110 ,3.906005 ,5.388615 ,
8.721716 ,8.875357 ,5.902633 ,5.821566 ,0.000000 ,7.584519 ,4.980863 ,6.500689 ,2.833213 ,7.059188 ,7.352377 ,2.975530 ,
3.837299 ,4.761319 ,1.131402 ,4.551769 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,0.000000 ,4.625953 ,4.151040 ,6.586723 ,
7.493985 ,6.640268 ,4.884316 ,0.000000 ,4.519612 ,4.298645 ,0.000000 ,7.911287 ,5.994709 ,4.480740 ,7.201618 ,5.085124 ,
5.917280 ,8.171006 ,7.923928 ,8.885040 ,7.063904 ,6.515009 ,4.117410 ,6.250168 ,6.563573 ,8.478120 ,0.000000 ,4.269697 ,
2.525729 ,0.000000 ,4.670958 ,3.572346 ,4.229749 ,5.709102 ,0.000000 ,4.105944 ,0.000000 ,5.217649 ,6.376387 ,2.674149 ,
5.093750 ,8.110007 ,7.615988 ,7.519041 ,5.631212 ,7.375632 ,4.990433 ,3.912023 ,6.855409 ,5.408964 ,5.607639 ,6.120737 ,
7.185008 ,7.449905 ,8.540851 ,8.529398 ,8.466152 ,7.608077 ,7.986301 ,2.980619 ,9.095535 ,8.420903 ,8.025026 ,5.921310 ,
5.986452 ,6.422110 ,5.830121 ,5.498806 ,6.689226 ,5.053695 ,6.601773 ,7.018133 ,7.742662 ,5.861071 ,7.381440 ,6.577444 ,
8.495132 ,8.823589 ,9.126415 ,7.851817 ,0.000000 ,7.340900 ,0.000000 ,0.000000 ,0.000000 ,8.711295 ,9.350041 ,0.000000 ,
0.000000 ,0.000000 ,5.776723 ,0.000000 ,0.000000 ,6.040493 ,5.625821 ,0.000000 ,6.680855 ,4.984291 ,0.000000 ,0.000000 ,
3.339322 ,5.110782 ,5.923721 ,6.903044 ,0.000000 ,4.291828 ,3.735286 ,9.188534 ,4.479607 ,6.339124 ,8.273872 ,1.722767 ,
7.236987 ,7.064674 ,8.598257 ,9.138984 ,10.463900 ,4.964940 ,6.652347 ,9.137243 ,0.000000 ,5.446306 ,8.385854 ,2.965273 ,
6.739337 ,7.443605 ,6.584101 ,5.690022 ,0.000000 ,8.037963 ,0.000000 ,5.498397 ,7.052288 ,6.843963 ,0.000000 ,5.605067 ,
6.000176 ,5.009301 ,0.000000 ,0.000000 ,6.292310 ,8.232998 ,0.000000 ,0.000000 ,8.035311 ,7.624228 ,8.464130 ,6.576609 ,
4.609162 ,7.243799 ,7.044295 ,8.031158 ,6.219198 ,3.054001 ,5.782286 ,6.924022 ,6.473119 ,6.043345 ,6.542040 ,8.781801 ,
8.095904 ,8.405546 ,7.392955 ,5.246498 ,6.520768 ,9.861024 ,5.974573 ,0.000000 ,6.297478 ,6.423734 ,5.451896 ,4.477337 ,
6.164157 ,7.937446 ,8.175182 ,7.071488 ,0.000000 ,0.000000 ,5.937800 ,0.000000 ,5.705115 ,0.000000 ,3.931826 ,3.139833 ,
0.000000 ,0.095310 ,4.514151 ,0.000000 ,7.684646 ,5.027165 ,4.643429 ,0.000000 ,6.200509 ,0.000000 ,4.582925 ,8.977475 ,
3.417727 ,4.222445 ,4.941642 ,8.638455 ,7.681653 ,7.306599 ,8.276191 ,8.569482 ,0.741937 ,4.235555 ,4.226834 ,5.716699 ,
6.396930 ,7.239215 ,6.240081 ,7.767730 ,7.922986 ,6.380631 ,0.000000 ,5.105945 ,7.726477 ,4.053523 ,2.624669 ,5.253843 ,
8.273923 ,5.070161 ,5.508578 ,6.545350 ,5.038250 ,0.000000 ,4.875960 ,6.201118 ,0.530628 ,0.000000 ,5.777962 ,8.164482 ,
6.302253 ,5.496348 ,8.042828 ,6.124246 ,7.394739 ,2.687847 ,5.342813 ,8.524268 ,8.204891 ,6.187031 ,5.860786 ,6.048317 ,
5.991715 ,7.213989 ,5.479805 ,0.000000 ,6.907655 ,6.926969 ,8.107027 ,7.313354 ,6.579112 ,4.568506 ,5.101085 ,5.195177 ,
3.648057 ,0.000000 ,5.007965 ,8.620832 ,5.212759 ,4.175925 ,2.772589 ,0.000000 ,4.890349 ,2.867899 ,5.023881 ,6.004134 ,
5.797880 ,5.673667 ,6.352978 ,5.204007 ,4.388257 ,5.319590 ,0.000000 ,6.869638 ,7.625693 ,5.324959 ,5.533389 ,5.259057 ,
5.673667 ,0.000000 ,7.304381 ,0.000000 ,6.135781 ,3.502550 ,0.000000 ,0.000000 ,7.279526 ,5.936480 ,0.000000 ,0.000000 ,
5.719984 ,7.112327 ,6.132096 ,4.733563 ,0.000000 ,4.799091 ,0.000000 ,7.721216 ,4.177459 ,2.533697 ,0.000000 ,0.000000 ,
0.000000 ,5.544787 ,7.451764 ,6.550509 ,6.593455 ,6.156767 ,4.722064 ,0.000000 ,6.000176 ,7.365876 ,0.000000 ,3.627004 ,
3.841601 ,5.258536 ,7.019118 ,7.610704 ,6.785588 ,0.000000 ,0.000000 ,2.557227 ,6.721787 ,4.935912 ,7.321850 ,6.283574 ,
5.544396 ,4.923624 ,5.940960 ,6.806940 ,0.993252 ,3.496508 ,7.478226 ,0.000000 ,4.970508 ,4.900076 ,0.000000 ,5.043425 ,
0.000000 ,6.254406 ,7.795688 ,7.591054 ,0.000000 ,5.188503 ,0.000000 ,6.283761 ,4.635699 ,0.000000 ,5.249652 ,0.000000 ,
0.000000 ,6.841402 ,0.000000 ,5.950122 ,0.000000 ,4.455509 ,7.377259 ,6.976721 ,0.000000 ,0.000000 ,5.598792 ,0.000000 ,
0.000000 ,0.000000 ,4.749271 ,4.442651 ,6.839369 ,7.012746 ,6.490268 ,7.918374 ,4.599152 ,6.595234 ,4.825109 ,5.236974 ,
0.000000 ,7.292814 ,4.494239 ,5.974318 ,7.958332 ,6.440149 ,3.673766 ,5.939908 ,0.000000 ,0.000000 ,8.835661 ,7.991558 ,
7.628761 ,5.697764 ,5.487283 ,8.634478 ,4.410371 ,8.000115 ,6.797494 ,6.455827 ,4.007333 ,4.623992 ,5.190175 ,9.976603 ,
6.213407 ,6.200103 ,4.465908 ,5.612763 ,5.894954 ,4.936630 ,6.055143 ,5.709102 ,5.080783 ,6.217604 ,4.539030 ,0.000000 ,
5.883044 ,6.561031 ,6.002899 ,5.775172 ,5.811141 ,6.048317 ,8.124062 ,0.000000 ,3.841601 ,0.000000 ,5.919969 ,6.557488 ,
7.243727 ,0.000000 ,3.072693 ,3.397858 ,6.406385 ,0.000000 ,0.000000 ,7.753108 ,6.777874 ,6.117657 ,7.068002 ,0.000000 ,
4.418841 ,4.280824 ,0.000000 ,4.588024 ,1.824549 ,0.000000 ,5.028475 ,3.068053 ,6.426650 ,0.000000 ,0.000000 ,7.122060 ,
7.561902 ,0.000000 ,6.318788 ,7.696894 ,0.000000 ,5.249652 ,4.055257 ,7.460950 ,0.000000 ,3.864931 ,4.918520 ,5.271973 ,
0.000000 ,3.478158 ,0.000000 ,0.000000 ,0.000000 ,4.816241 ,7.544808 ,7.304987 ,0.000000 ,5.732370 ,4.940213 ,4.511958 ,
4.725616 ,5.510602 ,1.648659 ,5.308763 ,8.238801 ,6.589064 ,7.172808 ,7.058328 ,7.717040 ,5.073297 ,2.509599 ,5.892749 ,
8.244465 ,7.139660 ,6.952729 ,6.057954 ,0.000000 ,9.540493 ,7.079184 ,7.521589 ,7.653922 ,6.382493 ,9.095233 ,7.244799 ,
8.046677 ,7.153286 ,9.186181 ,5.284218 ,8.748495 ,0.000000 ,0.000000 ,5.247550 ,8.275300 ,6.929517 ,6.507725 ,5.919163 ,
6.817721 ,5.597681 ,6.606920 ,5.401325 ,5.775793 ,5.632644 ,4.781641 ,6.908655 ,5.424509 ,4.867534 ,6.473119 ,7.571474 ,
5.476464 ,4.716712 ,4.748404 ,6.789310 ,6.722148 ,4.731803 ,5.570632 ,6.823068 ,7.071488 ,4.575741 ,6.953780 ,5.046002 ,
6.296188 ,8.004766 ,6.463341 ,3.873282 ,4.354141 ,0.000000 ,3.779634 ,0.000000 ,6.234018 ,5.224671 ,4.330733 ,5.493061 ,
5.583120 ,6.871299 ,5.554122 ,7.918228 ,0.000000 ,4.334673 ,5.571774 ,5.843255};







	TrainModelEM(&trained_hmm, cancer_exp_vec, 1104,  max_iters, num_starts, tolerance, best_model_score);
	PrintModel(&trained_hmm);


	// free all allocated buffers
	delete x_vec;
	delete y_vec;
	delete y_vec_int;
	delete trained_x_vec;
	delete loc_vec;
	delete scale;
	delete scale_cum;
	delete scale_exp;

	for(i = 0; i < my_hmm.x_dim; i++)
	{
		delete y_cond_tabs[i];
		delete alpha[i];
		delete beta[i];
		delete gamma[i];

		for(j = 0; j < my_hmm.y_dim; j++)
		{
			delete gammagamma[i][j];
			delete count_occur[i][j];
			delete y_mu_square_exp_vecs[i][j];
		}
	}


	// problemmm !!! maybe need a different structure for the data ??? 
	delete my_data.y_vec;
	delete my_data.y_vec_int;  
	delete my_data.x_vec;
	delete my_data.mix_vec;
	delete my_data.loc_vec;
	delete my_data.loc_diff_vec;

#endif

	return 0; 
}





