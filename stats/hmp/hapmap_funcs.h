#ifndef _HAPMAP_FUNCS_
#define _HAPMAP_FUNCS_


#define NUM_TRIOS 90 

long ReadGenotypeFile(FILE *genotype_f, char *genotype_p,  // input
					  snp_hapmap_data *snp_hapmap);
long CountFileLines(FILE *f, char *p);
long PopCount(word w);

#endif