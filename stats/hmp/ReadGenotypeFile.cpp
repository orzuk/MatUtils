#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "general.h"
#include "hmm_chrom_funcs.h"
#include "hapmap_funcs.h"


///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For HAPMAP Chromosome project - 11.09.2006
///////////////////////////////////////////////////////////////////////////////////////
//#undef MEX_TO_MATLAB
//#ifdef MEX_TO_MATLAB
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  
	int nRows, nCols, i, j, K, t, max_iters, num_starts;

	long seq_len, L; 
	long x_dim, y_dim, place_flag, miss_data;
	double tolerance;
  
	double *in, *out;

	
	FILE *genotype_f;  
//	char genotype_p[150] = "E:\\Research\\HMM_Chromosome\\SNP\\HAPMAP\\Genotype_data\\genotypes_chr"; // 22_CEU_r21_nr_fwd.txt"; // later we must use the chromosome as input here!!!!!  
////	char genotype_p[150] = "genotypes_chr"; // 22_CEU_r21_nr_fwd.txt"; // later we must use the chromosome as input here!!!!!  
	char *genotype_p;

	 /* Check for proper number of arguments. Should be One (1) */
	if(nrhs != 1/*3*/) 
		mexErrMsgTxt("Usage: [dat] [K]");

  /* Parameters should be : 
	1. Chromosome number
	2. Population String
	3. Hapmap Version string

	2. Offsprings Indexes(???) - We Ignore these ones when doing the statistics - CANCELLED !!! 
	*/


/*** New: We cut the bullshit and put it in the Matlab. Here, all that we get is the file name.
Parameters should be : 
	1. File name string
	***/



  // Get the directory we work on ? NO! We need to go to the database directory in Matlab ourselves !!!  
//	printf("Start All This Shit  \n");


/***	

	in = mxGetPr(prhs[0]);  // Get the chromosome number 
	long chrom = long(in[0]);  	
	char *chrom_p = new char[2];
	itoa(chrom, chrom_p, 10);
//	printf("The chrom is : %s\n", chrom_p);

	if(chrom == 23) // The X chromosome
		chrom_p = "X";
	if(chrom == 24) // The Y chromosome
		chrom_p = "Y";
	if(chrom == 25) // The Mitochondria
		chrom_p = "M";
//	printf("Now The chrom is : %s\n", chrom_p);


//	in = mxArrayToString(prhs[1]);  // Get the population  string 

	char *population_p; // = new char[20];
	population_p = mxArrayToString(prhs[1]);
	char *hap_version_p; // = new char[20];
	hap_version_p = mxArrayToString(prhs[2]);

	printf("Population Is %s  \n", population_p);


	// now try to 'paste' 
	strcat(genotype_p, chrom_p); 
	strcat(genotype_p, "_");
	strcat(genotype_p, population_p);
	strcat(genotype_p, "_"); // previousely it was "_r" - removed 
	strcat(genotype_p, hap_version_p);
	strcat(genotype_p, "_nr_fwd.txt");
//	CEU_r21_nr_fwd.txt"); 
		printf("File name nnnow: %s\n", genotype_p);
	delete chrom_p;
***/
	genotype_p = mxArrayToString(prhs[0]); // Get the file name


	snp_hapmap_data snp_hapmap;

//	printf("Just before counting lines ...\n");
//	printf("File name nnnow: %s\n", genotype_p);
	snp_hapmap.num_lines = CountFileLines(genotype_f, genotype_p);
	
//	printf("File contains %d lines\n", snp_hapmap.num_lines);

	
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

//	printf("Done reading file\n");

	//	printf("num lines outside is %d\n", snp_hapmap.num_lines);



	// Now start copying the output matrices  : 
	plhs[0] = mxCreateDoubleMatrix(6, snp_hapmap.num_lines, mxREAL);  // The SNPs data 
	out = mxGetPr(plhs[0]); // output matrix
	for(i=0; i<snp_hapmap.num_lines; i++) 
		for(j=0; j<6; j++) 
			out[6*i+j] = snp_hapmap.data[i][j];

	plhs[1] = mxCreateDoubleMatrix(snp_hapmap.num_lines, 1, mxREAL);  // The SNPs bases 
	out = mxGetPr(plhs[1]); // output vector
	for(i = 0;i < snp_hapmap.num_lines; i++)
			out[i] = snp_hapmap.bases[i];

	plhs[2] = mxCreateDoubleMatrix(snp_hapmap.num_lines, 1, mxREAL);  // The SNPs frequencies
	out = mxGetPr(plhs[2]); // output vector
	for(i = 0;i < snp_hapmap.num_lines; i++)
			out[i] = snp_hapmap.freqs[i];

	plhs[3] = mxCreateDoubleMatrix(snp_hapmap.num_lines, 1, mxREAL);  // The SNPs hetrozigocity frequencies
	out = mxGetPr(plhs[3]); // output vector
	for(i = 0;i < snp_hapmap.num_lines; i++)
			out[i] = snp_hapmap.hetros[i];

	plhs[4] = mxCreateDoubleMatrix(6, snp_hapmap.num_lines, mxREAL);  // The SNPs bad calls vec 
	out = mxGetPr(plhs[4]); // output matrix
	for(i=0; i<snp_hapmap.num_lines; i++) 
		for(j=0; j<6; j++) 
			out[6*i+j] = snp_hapmap.bad_calls[i][j];

	plhs[5] = mxCreateDoubleMatrix(snp_hapmap.num_lines, 1, mxREAL);  // The SNPs chromosomal locations
	out = mxGetPr(plhs[5]); // output vector
	for(i = 0;i < snp_hapmap.num_lines; i++)
			out[i] = snp_hapmap.chrom_locs[i];


	// Copy the name strings	
	plhs[6] = mxCreateCharMatrixFromStrings(snp_hapmap.num_lines, (const char **)(snp_hapmap.names));


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

}
//#endif // MEX_TO_MATLAB



