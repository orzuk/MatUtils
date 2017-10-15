#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "general.h"
#include "hapmap_funcs.h"


// Functions specially dealing with hapmap data, SNP chips etc.




/////////////////////////////////////////////////////////////////////
// Reading a genotype file and extracting all the relevant SNPs
// Output C = A*B is r*n. We assume C is allocated 
// we assume the maximal dimentions are MAX_X_VALS 
//
//
long ReadGenotypeFile(FILE *genotype_f, char *genotype_p,  // input
					  snp_hapmap_data *snp_hapmap)				
{
	long i, j, s; //, num_lines;
//	char c;	
	char tmp_str[1000], tmp_short_str[50];
//	int first_space;
	char spc = ' ';
	char QC_plus = 'QC+';
	char *first_space_p, *second_space_p, *third_space_p, *fourth_space_p; 
	char *start_seqs_p; // *fff_p,
	long first_space_num;
	long bad_calls_num, bad_calls_hetro_num;
//	long location;
	

	return 9999;

	// Open again and start reading
	genotype_f = fopen(genotype_p, "r");
	for(i = 0; i < snp_hapmap->num_lines; i++)	
	{	
		snp_hapmap->freqs[i] = 0;
		fgets(tmp_str, 1000, genotype_f); 

		// Get the SNPs name
		first_space_p = strchr(tmp_str, spc ); 
		second_space_p = strchr(first_space_p+1, spc ); 
		third_space_p = strchr(second_space_p+1, spc ); 
		fourth_space_p = strchr(third_space_p+1, spc ); 
		first_space_num = first_space_p-tmp_str+1;

//		if(i == 64)
//			printf("Line Is: %s\n", tmp_str);

		for(j=0; j < 6; j++)
		{
			snp_hapmap->data[i][j] = 0UL; // Make everything zero (the Major allele)
			snp_hapmap->bad_calls[i][j] = 0UL; // Make bad-calls zero (assume good calls)
		}

		strncpy (snp_hapmap->names[i], tmp_str, first_space_num-1); 
		snp_hapmap->names[i][first_space_num-1] = NULL;

		
		// Get the SNPs chromosomal location
		start_seqs_p = strstr( tmp_str, "chr");
		if(start_seqs_p != NULL)		
		{
			strncpy (tmp_short_str, start_seqs_p+5, fourth_space_p-(start_seqs_p+5)); 
			tmp_short_str[fourth_space_p-(start_seqs_p+5)+1] = NULL;
			snp_hapmap->chrom_locs[i] = atoi(tmp_short_str);
		}


		// Get the SNPs Letters
		//printf("Letter1 is %c\n", first_space_p[1]);
		switch (first_space_p[1])
		{			
			case 'A' : 	snp_hapmap->bases[i] = 0; break;
			case 'T' :	snp_hapmap->bases[i] = 1; break;
			case 'C' :  snp_hapmap->bases[i] = 2; break;
			case 'G' :  snp_hapmap->bases[i] = 3; break;
			default : snp_hapmap->bases[i] = 99; break;
		}

		//printf("Letter3 is %c\n", first_space_p[3]);
		switch (first_space_p[3])
		{

			case 'A' : snp_hapmap->bases[i] += (4*0); break;
			case 'T' : snp_hapmap->bases[i] += (4*1); break;
			case 'C' : snp_hapmap->bases[i] += (4*2); break;
			case 'G' : snp_hapmap->bases[i] += (4*3); break;
			default : snp_hapmap->bases[i] = 99; break;
		}

		// Now get the data		
		start_seqs_p = strstr( tmp_str, "QC+");
		if(start_seqs_p != NULL)
		{
			start_seqs_p = start_seqs_p+4; 
			for(s=0; s<NUM_TRIOS; s++)  // loop over 90 patients ... 
			{
				// Set one only if we're in the minor allele
				if(start_seqs_p[3*s] == first_space_p[3])
					snp_hapmap->data[i][s/30] ^= (1UL << (s%30));
				if(start_seqs_p[3*s+1] == first_space_p[3])
					snp_hapmap->data[i][(s/30)+3] ^= (1UL << (s%30));

				// Find the bad calls
				if(start_seqs_p[3*s] == 'N')
					snp_hapmap->bad_calls[i][s/30] ^= (1UL << (s%30));
				if(start_seqs_p[3*s+1] == 'N')
					snp_hapmap->bad_calls[i][(s/30)+3] ^= (1UL << (s%30));
			}
		}

		// Now calculate the frequency
		bad_calls_num = 0;
		for(j=0;j<6;j++)
		{
			bad_calls_num += PopCount(snp_hapmap->bad_calls[i][j]);
			snp_hapmap->freqs[i] += PopCount(snp_hapmap->data[i][j] & (~snp_hapmap->bad_calls[i][j]));
			if(i == 63)
				printf("Bads are %d\n", bad_calls_num);
		}
		bad_calls_hetro_num = 0;
		for(j=0; j<3;j++)
		{
			bad_calls_hetro_num += PopCount(snp_hapmap->bad_calls[i][j] | snp_hapmap->bad_calls[i][j+3]);
			snp_hapmap->hetros[i] += PopCount( (snp_hapmap->data[i][j] ^ snp_hapmap->data[i][j+3]) & 
												(~(snp_hapmap->bad_calls[i][j] | snp_hapmap->bad_calls[i][j+3])) ); // Take only when two alleles are different
		}
		snp_hapmap->freqs[i] = 1.0-snp_hapmap->freqs[i]/(180.0-1*bad_calls_num); // take the Major frequency 
		snp_hapmap->hetros[i] /= (90.0-bad_calls_hetro_num); // take the Major frequency 

		// Now actually 'READ' the line and extract the meaningfull information out of it	
		if(start_seqs_p != NULL)
		if(i == 63)
		{
			printf("Line is:  %s", tmp_str);
			printf("First space is on %d\n", first_space_num);
			printf("SNP name : %s\n", snp_hapmap->names[i]);
			printf("Letters are: %c and %c\n", first_space_p[1], first_space_p[3]);  
			printf("First data letters %c%c  %c%c\n", start_seqs_p[0],start_seqs_p[1],start_seqs_p[3*1],start_seqs_p[3*1+1]);
			printf("Freq is: %lf\n", snp_hapmap->freqs[i]);			
			printf("Data is:  %s\n\n", start_seqs_p);
		}
	}
	return snp_hapmap->num_lines;
}



// Simply the count the number of lines in a file
long CountFileLines(FILE *f, char *p)
{
	long i, num_lines=999;
	char *fff_p;
	char tmp_str[1000]; // assume this is the maximal line length - is this correct ??? 

	// open the file to enable reading from it
	f = fopen(p, "r");

	// find length of file	
	i=0;
	while(fff_p != NULL)  // should be in ... 
/////	while(i == 1)
	{
		fff_p = fgets(tmp_str, 1000, f);
		i++;
/////		if(i==1000) // stop early .. (WRONG !!!)
////		{
//////			fclose(f);
/////			return num_lines;
/////		}

//		if((i%1000)==0)
//			printf("I is %d\n", i);
	}

	num_lines = i-1;

	fclose(f);

	return num_lines;
}




// Count the number of on bits in a word
long PopCount(word w)
{
	long pop=0;
	long i;

	for(i=0;i<8;i++)
		pop+=pop_tab[(w >> (4*i))&0xF];
	return pop;


}

