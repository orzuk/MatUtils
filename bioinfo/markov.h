// File markov.h - Macros, prototypes and constant of the DNA package
#ifndef _MARKOV_H_
#define _MARKOV_H_

long transfer_weights(double w[4][4][MAX_L], long m[4][4][MAX_L], double res, long l);
double calc_thresh_prob(long m[4][4][MAX_L], long thresh, long l, double *probs, double *dist, long *max_val, long *min_val);
double calc_thresh_prob_from_dist(long thresh,  double *dist, long strands, long max_val, long l, long trials);
double my_power(double m, long d);

#endif
