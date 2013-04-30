#ifndef _DNA_UTILS_H_
#define _DNA_UTILS_H_


void print_gene(gene g);
void print_dna(char *seq, long seq_len);
void print_bin(word *a, long L);
void print_quad(word *a, long L);
void print_rna(word *a, long L);
//void print_packed_dna(word *seqs, long len);
void print_packed_dna(word *seqs, long len, long print_mode);

long find_gene_length(FILE *promoter_f, char *promoter_p);
long read_gene(FILE *promoter_f, char *promoter_p, gene *g);
long read_all_genes(FILE *promoter_f, char *promoter_p, gene *g);
long read_expression(FILE *expression_f, char *expression_p, gene *g, char labels[MAX_SAMPLES][MAX_NAME_LENGTH] );
long calc_correlation(gene *g, double *expression);
void calc_reverse_complement(char *seq, long seq_len, char *comp_seq, long *comp_seq_len);
long expand_seq(char *seq, long seq_len , int *seq_exp, long *seq_exp_len);


// We currently assume a markov model of size up to 6
void markov_model_collect_statistics(char *seq, long seq_len,  long counts[4][4][4][4][4][4]);
void print_probs(double *probs, long dim);
void count_to_probs(long counts[4][4][4][4][4][4], double probs[4][4][4][4][4][4]);
long update_count(long *counts, char *seq, long seq_len);
void counts_to_probs(long *counts, long total_count, long dim, double derich, double *probs);
void generate_random_seq(char *seq, long seq_len, double *probs, long dim);
void generate_random_seq_word(word *seq, long seq_len, double *probs, long dim);
void probs_to_weights(double *probs, long dim, double *weights);
long sort_all_genes(gene *g);
void probs_to_markov_probs(double *probs, double *markov_probs, long dim);
void read_markov_probs(FILE *probs_f, char *probs_p, double *probs);
void read_threshes(FILE *threshs_f, char *threshs_p, double genes_threshes[L_TF][L_GENES]);
void read_TF_matrices(FILE *tfname_f, char *tfname_p, FILE *gename_f, char *gename_p, 
					  long pssm[L_TF][4][MAX_L], long pssm_len[MAX_L], 
					  char TF_NAMES[L_TF][MAX_NAME_LENGTH], char gene_names[L_GENES][MAX_NAME_LENGTH]);
void pssm_and_background_to_weights(long pssm[4][MAX_L], long pssm_len, double derich, double *back_probs, 
									double w[4][4][MAX_L]);
double calc_site_score(double w[4][4][MAX_L], long l, long seq[MAX_L]);
void write_probs(FILE *probs_f, char *probs_p, double probs[L_TF][L_GENES]);
void read_lengths(FILE *lengths_f, char *lengths_p, long genes_lengths[L_GENES]);
long calc_all_sites_scores(double weights[4][MAX_L], long L, word *seqs[MAX_NUM_SEQS], 
						   double *scores[MAX_NUM_SEQS], long num_seqs, long *seqs_lens); 
long calc_all_sites_scores_single(float weights[4][MAX_L], long L, word *seqs[MAX_NUM_SEQS], 
						   float *scores[MAX_NUM_SEQS], long num_seqs, long *seqs_lens);
long get_scores_above_threshold(long is_double, double *scores[MAX_NUM_SEQS], float *scores_single[MAX_NUM_SEQS], 
								long num_seqs, long *seqs_lens, double threshold, float threshold_single,
								long remove_overlap_distance, long smart_thresh, 
								long *out_x, long *out_y, 
								double *out_scores, float *out_scores_single); // Those are allocated OUTSIDE the function
long match_closest_vals(double *closest_ind, double *closest_dist, double *a, double *b, long a_len, long b_len);
long intervals_intersect(double *intersect_start_pos, double *intersect_end_pos, double *intersect_inds1, double *intersect_inds2, long *inter_ctr,
		double *start_pos1, double *end_pos1, double *start_pos2, double *end_pos2, long n1, long n2, long is_sorted);
long extract_sub_kmers(long L, word *seqs[MAX_NUM_SEQS], long *seqs_lens, long num_seqs, long unique_flag, 
					   word *kmers[MAX_L], long *kmer_kmer_inds, long *kmer_seq_inds, long *kmer_position_inds, 					   
					   long *num_unique_kmers, long *num_sparse_matrix_elements, long hash_flag, 
						long input_position_flag, long num_coordinates);
long add_noise_to_kmers(long L, long num_kmers, word *kmers[MAX_L], double noise_table[MAX_L][4][4], 
						word *noisy_kmers[MAX_L]);


#define MEX_TO_MATLAB
/*
#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[]);
#endif // MEX_TO_MATLAB
*/

#endif
