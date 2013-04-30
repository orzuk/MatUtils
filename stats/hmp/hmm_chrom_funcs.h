#ifndef _HMM_CHROM_FUNCS_
#define _HMM_CHROM_FUNCS_

long split(double *vals, long first, long last, long *indexes);
long quicksort(double *vals, long first, long last, long *indexes);
long DoQuicksort(double *vals, long len, long *indexes);
long DoOrder(double *vals, long len, long *indexes);

long MatrixMultiply(double A[MAX_X_VALS][MAX_X_VALS], long m, long n, double B[MAX_X_VALS][MAX_X_VALS], long r, 
					double C[MAX_X_VALS][MAX_X_VALS]);
long MatrixStationaryVec(double A[MAX_X_VALS][MAX_X_VALS], long n, double pi[MAX_X_VALS]);
long MatrixTranspose(double A[MAX_X_VALS][MAX_X_VALS], long m, long n, 
					 double A_t[MAX_X_VALS][MAX_X_VALS]);
long GSolve(double a[MAX_X_VALS][MAX_X_VALS],long n,double x[MAX_X_VALS]);
					

double MixtureOfGaussiansGivenInit(double *x, long num_points, long num_of_Gaussians, long num_of_iterations, 
		 double *init_p, double *init_m, double *init_s, 
		 double *p, double *m, double *s, double *log_like);
float MixtureOfGaussiansGivenInitSingle(float *x, long num_points, long num_of_Gaussians, long num_of_iterations, 
		 float *init_p, float *init_m, float *init_s, 
		 float *p, float *m, float *s, float *log_like);



double randr();
long PrintVec(long *vec, long len);
long PrintDoubleVec(double *vec, long len);
long Geometric(double p);
double Gaussian(double mu,double sigma);
double ComputeJointHMMLogLikelihood(hmm_model *hmm, hmm_data *data);
long Viterbi( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_X_VALS],  double *y_cond_tabsB[MAX_X_VALS],
			 long *opt_x_vec_out, double *x_loglike, double *x_place_loglike, double *y_loglike);
long InitilizeData( hmm_data *data, long seq_len, long y_type, long miss_data, double dont_miss_prob, long copy_data, 
				   double *y_vec, long *y_vec_int, long *loc_vec);
long InitilizeModel(hmm_model *hmm, long x_dim, long y_dim, long y_type, long use_bounds, long place_flag, 
					long special_models_flag, long miss_data, long fold_change_flag);
long ComputeModelAuxillaryParameters(hmm_model *hmm);
long CompareTwoModels(hmm_model *hmm1, hmm_model *hmm2);
long SimulateSequenceFromModel(hmm_model *hmm, hmm_data *data);
long PrintModel(hmm_model *hmm);
long forward( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_X_VALS], 
			  double *alpha[MAX_X_VALS], double *scale, double *scale_cum,  double *scale_exp, double *y_vec_log_prob_out);
long backward( hmm_model *hmm, hmm_data *data, double *scale_exp, double *y_cond_tabs[MAX_X_VALS],
			  double *beta[MAX_X_VALS]); 
long ComputeMarginalGamma(hmm_model *hmm, hmm_data *data, double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS],
						 double *alpha[MAX_X_VALS], double *beta[MAX_X_VALS], double *scale_exp, double *y_cond_tabs[MAX_X_VALS],
						 double *gamma[MAX_X_VALS], double *phi[MAX_X_VALS][MAX_Y_VALS]);
long ComputePsi(hmm_model *hmm, hmm_data *data, double *alpha[MAX_X_VALS], double *beta[MAX_X_VALS], double *scale, double *y_cond_tabs[MAX_X_VALS],
				 double *psi[MAX_X_VALS][MAX_X_VALS]);
long Compute_y_cond_tabs(hmm_model *hmm,  hmm_data *data, 
					double *y_cond_tabs[MAX_X_VALS], 
					double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS]);
long UpdateModelParams(hmm_model *hmm, hmm_data *data, double *y_mu_square_exp_vecs[MAX_X_VALS][MAX_Y_VALS],
						 double *gamma[MAX_X_VALS], double *phi[MAX_X_VALS][MAX_Y_VALS],
						 double *psi[MAX_X_VALS][MAX_X_VALS]);
long SmartStartParams(hmm_model *hmm, double *y_vec, long seq_len);
long TrainModelEM(hmm_model *hmm, hmm_data *data, long max_iters, long num_starts, double tolerance, 
				  double *best_model_score);
long CopyModel(hmm_model *hmm1, hmm_model *hmm2, hmm_data *data); 
long PermuteModelIncreasingMean(hmm_model *hmm, hmm_data *data);
long ComputeKLDistance(hmm_model *hmm1, hmm_model *hmm2, hmm_data *data, long iters, double *KL_dist);

// special functions for the 'complicated' SNPs case
//double ComputeJointHMMLogLikelihoodSNPs(hmm_model *hmm, hmm_data *data);
double ComputeJointHMMLogLikelihoodSNPs(hmm_model *hmm, hmm_data *data,
										double *x_loglike, double *place_x_loglike, double *y_loglike);
double TrainModelEMSNPs(hmm_model *hmm, hmm_data *data, long max_iters, long num_starts, double tolerance, 
				  double *best_model_score, double *tmp_vec);
long forwardSNPs( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS], 
			  double *alpha[36], double *scale, double *scale_cum,  double *scale_exp, double *y_vec_log_prob_out);
long backwardSNPs( hmm_model *hmm, hmm_data *data, double *scale_exp, 
				   double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS],
				   double *beta[36]);
long ViterbiSNPs( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_2D_X_VALS],  double *y_cond_tabsB[MAX_2D_X_VALS],
			  long *opt_x_vec_out); // , double *x_loglike, double *x_place_loglike, double *y_loglike);
long ComputeMarginalGammaSNPs(hmm_model *hmm,  hmm_data *data, 
						 double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS], double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS],
						 double *alpha[36], double *beta[36], double *scale_exp, 
						 double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS],
						 double *gamma[36], double *phi[36][MAX_Y_VALS]);
double UpdateModelParamsSNPs(hmm_model *hmm, hmm_data *data, 
						   double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS],
						   double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS],
						   double *gamma[36], double *phi[36][MAX_Y_VALS],
						   double psi_sum[36][36]);
long ComputePsiSNPs(hmm_model *hmm, hmm_data *data, double *alpha[36], double *beta[36], 
				double *scale, double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS],
				double psi_sum[36][36]);

long Compute_y_cond_tabsSNPs(hmm_model *hmm,  hmm_data *data, 
					double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS], 
					double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS], 
					double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS]);


long ComputeKLDistanceSNPs(hmm_model *hmm1, hmm_model *hmm2, hmm_data *data, long iters, double *KL_dist);


long NEWCompute_y_cond_tabsSNPs(hmm_model *hmm,  hmm_data *data, 
					double *y_cond_tabs[36], double *y_cond_tabsB[36], 
					double *y_mu_square_exp_vecs[MAX_2D_X_VALS][MAX_Y_VALS], 
					double *y_mu_square_exp_vecsB[MAX_2D_X_VALS][MAX_Y_VALS]);

long SLOWforwardSNPs( hmm_model *hmm, hmm_data *data, double *y_cond_tabs[MAX_2D_X_VALS], double *y_cond_tabsB[MAX_2D_X_VALS], 
			  double *alpha[36], double *scale, double *scale_cum,  double *scale_exp, double *y_vec_log_prob_out);

#endif