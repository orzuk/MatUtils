#include <stdlib.h>
#include <stdio.h>
// #include <iostream.h>
#include <time.h>
#include <math.h>
#include <string.h>

//#include "mex.h"
#include "general.h"
#include "dna_utils_new.h"
#include "markov.h"





/////////////////
// A list of functions dealing with DNA
// Note : Here we deal with DNA in COMPRESSED form. 
// This will be the standard from now on !!!!!!!
// We will use chars of 8 bits. The 2 lsb will only be considered
// Each 2 bits will be one DNA base. 
// This will be : (from Now on !!!!)
// A - 00
// T - 01
// C - 10
// G - 11
// a DNA sequence will always be an array (of chars ???), with the first char saying the length !!!
// There should always be enough place allocated
// This causes a 4-times extra memory, but what the f***
/////////////////


// Print a gene gsdafsd
void print_gene(gene g)
{
	long i;

	printf("Locus = %d Contig = %d gene name = %s   gene symbol = %s afyname = %s \n",   
		g.locus, g.NT, g.name, g.symbol, g.afy_name);

	printf("Expression : ");
	for(i = 0; i < g.exp_len; i++)
		printf("%lf, ", g.exp_vec[i]);
	printf("\n");

	print_dna(g.seq, g.seq_len);

}



// Print the DNA sequence
void print_dna(char *seq, long seq_len)
{
	

	printf("seq len : %ld\n seq : ", seq_len);

	
/***	
	long i;
	for(i = 0; i < seq_len; i++)
	{
		if( (seq[i]&0x3) == 0)
			printf("A");
		else if( (seq[i]&0x3) == 1)
			printf("T");
		else if( (seq[i]&0x3) == 2)
			printf("C");
		else 
			printf("G");

		if((i%40) == 0 )
			printf("\n");
				

	}
***/
	printf("%s\n", seq);

}








////////////////////////////////////////////
//
//  function : find_gene_length
// 
//  input    : promoter_f, promoter_p
//  output   : g
//
//  purpose  : The function finds the length of 
//			   a gene from the input file
//			   
//
////////////////////////////////////////////
long find_gene_length(FILE *promoter_f, char *promoter_p)
{
	long i;

	char c;
	


	// open the file to enable reading from it
	promoter_f = fopen(promoter_p, "r");
	c = ' ';

	printf("Finding Length ....\n");

	// find length of file	
	i = 0;
	while(fscanf(promoter_f, "%c", &c) != EOF)	
		i++;
	

	fclose(promoter_f);

	return i;

}



////////////////////////////////////////////
//
//  function : read_gene
// 
//  input    : promoter_f, promoter_p
//  output   : g
//
//  purpose  : The function reads a gene from  
//			   a file.
//
////////////////////////////////////////////
long read_gene(FILE *promoter_f, char *promoter_p, gene *g)
{
	long i;

	char c;
	
	


	// open the file to enable reading from it
	promoter_f = fopen(promoter_p, "r");


//	printf("Start Reading The gene ...\n");

	c = ' ';

	printf("Start Reading The gene ... \n");
	while( c  == ' ')
		fscanf(promoter_f, "%c", &c);

	printf("Finished spaces Reading The gene ...\n");



	// Now reading the parameters (they must be in the correct order in the file)
	fscanf(promoter_f, "ocus=%d ", &g->locus);
	fscanf(promoter_f, "Contig=NT_%d", &g->NT);



	i = 0;
	c = ' ';
	while( int(c) != EOLN)
	{
		fscanf(promoter_f, "%c", &c);
		g->name[i++] = c;

	}
	g->name[i] = '\0';

	
	printf("Finished name fields Reading The gene ...\n");
	printf("Locus = %d Contig = %d gene name = %s\n",   g->locus, g->NT, g->name);
	
//	c == ' ';




	i = 0;
	while( fscanf(promoter_f, "%c", &c) != EOF )
	{
		g->seq[i++] = c;
		
	}
	g->seq[i] = '\0';
	printf("Gene was read successfully before closing\n"); 

	fclose(promoter_f);
	printf("Gene was read successfully\n"); 



	return 0;

}







////////////////////////////////////////////
//
//  function : read_all_genes
// 
//  input    : promoter_f, promoter_p
//  output   : g
//
//  purpose  : The function reads all genes from
//			   a file 
//			   
//
////////////////////////////////////////////
long read_all_genes(FILE *promoter_f, char *promoter_p, gene *g)
{
	long i;
	
	long genes_num;

	char c;
	
	


	// open the file to enable reading from it
	promoter_f = fopen(promoter_p, "r");


	printf("Start Reading Whole The genes ...\n");


	c = ' ';


	for(genes_num = 0; genes_num < MAX_GENES_NUM; genes_num++)
	{
	
		while(c  == ' ')
		{
			fscanf(promoter_f, "%c", &c);
			printf("c = %c \n", c);
		}




		// Now reading the parameters (they must be in the correct order in the file)
		fscanf(promoter_f, "Locus=%d ", &(g[genes_num].locus) );
		fscanf(promoter_f, "Contig=NT_%d", &(g[genes_num].NT) );



		i = 0;
		c = ' ';
		while( c == ' ')
			fscanf(promoter_f, "%c", &c);


		// check if we are in complement ...
		if(c == '*')
		{
			fscanf(promoter_f, "complement", &c);
			c = ' ';
			while( c == ' ')
				fscanf(promoter_f, "%c", &c);
		}



		while( (c != ' ') && (int(c) != EOLN) )
		{
			g[genes_num].name[i++] = c;
			fscanf(promoter_f, "%c", &c);
		}
		g[genes_num].name[i] = '\0';


		while( int(c) != EOLN)
			fscanf(promoter_f, "%c", &c);


	
		if( (genes_num%1000) == 0 ) 
			printf("Locus = %d Contig = %d name = %s ",   g[genes_num].locus, g[genes_num].NT, g[genes_num].name);
	
//		c == ' ';




		i = 1;
		g[genes_num].seq[0] = c; 
		while( g[genes_num].seq[i-1] /*c*/ != '>'  ) 		
		{
			fscanf(promoter_f, "%c", &(g[genes_num].seq[i++]));
//			fscanf(promoter_f, "%c", &c);
//			g[genes_num].seq[i++] = c;
			
		}
		g[genes_num].seq[i] = '\0';


		if( (genes_num%1000) == 0 ) 
			printf("Read g # %d \n", genes_num);



	}
	
	
	 fclose(promoter_f);


	 printf("\nAll Genes were read successfully\n"); 



	return 0;

}






////////////////////////////////////////////
//
//  function : read_expression
// 
//  input    : expression_f, expression_p
//  output   : g
//
//  purpose  : The function reads expression data
//			   from a file			   
//
/////////////////////////////////////////////
long read_expression(FILE *expression_f, char *expression_p, gene *g, char labels[MAX_SAMPLES][MAX_NAME_LENGTH] )
{
	long i, j, k;

	char c;
	
	


	// open the file to enable reading from it
	expression_f = fopen(expression_p, "r");


	printf("Start Reading The expression ...\n");

	c = 'k';

	printf("Start Reading The labels ... \n");

	i = 0;
	for(i = 0; i < RAMAS_SAMPLES; i++)
	{
		j = 0;

		while( int(c) == TAB) 
		{
			fscanf(expression_f, "%c", &c);

//			if(i >= RAMAS_SAMPLES-2)
//				printf("TREAD  %c %d\n", c, int(c));
		}

		while( (int(c) != TAB) && (int(c) != EOLN) ) // 9 is TAB	10 is EOLN

		{
	//		if( int(c) == EOLN )
	//			break;
			labels[i][j++] = c;
			
			fscanf(expression_f, "%c", &c);				

//			if(i >= RAMAS_SAMPLES-2)
//				printf("RREAD  %c %d\n", c, int(c));

		}
		labels[i][j] = '\0';
	}

	printf("Finished Reading The labels  ...\n");




//exit(99);

	c = 'k';
	for(k = 0; k < EXP_GENES_NUM; k++)
	{

		// Read the first two fields
		while( (int(c) == TAB) || (int(c) == EOLN) )
			fscanf(expression_f, "%c", &c);
	
		
		i = 0;
		while( int(c) != TAB) // 9 is TAB	
		{		
		//	g[genes_num].seq[i] = '\0';
			g[k].afy_name[i++] = c;
					
			fscanf(expression_f, "%c", &c);				
		}	
		g[k].afy_name[i] = '\0';


		while( int(c) == TAB) 
			fscanf(expression_f, "%c", &c);
	
		i = 0;
		while( int(c) != TAB) // 9 is TAB	
		{			
			g[k].symbol[i++] = c;
					
			fscanf(expression_f, "%c", &c);				
		}	
		g[k].symbol[i] = '\0';



		// Now read the numbers	
		for(i = 0; i < RAMAS_SAMPLES-2; i++)
			fscanf(expression_f, "%lf", &(g[k].exp_vec[i]));

	}


/*****************************************************************************************
   In Comment !!!
		{
			printf("i %d\n", i);
			j = 0;

			while( int(c) == TAB) 
				fscanf(expression_f, "%c", &c);
	
			while( int(c) != TAB) // 9 is TAB	

			{
				if( int(c) == EOLN )
					break;
				labels[i][j++] = c;
				
				fscanf(expression_f, "%c", &c);				
			}
			labels[i][j] = '\0';
		}	

	}

	while( int(c) != EOLN)
	{
		printf("i %d\n", i);
		j = 0;
		while( c  != ' ')
		{
	//		printf("scanning ...\n");
			fscanf(expression_f, "%c", &c);

			if( int(c) == EOLN )
				break;
	//		printf("READ %c\n", c);
			labels[i][j++] = c;
		}
		labels[i][j] = '\0';
		i++; 
	}

/***
	// Now reading the parameters (they must be in the correct order in the file)
	fscanf(promoter_f, "ocus=%d ", &g->locus);
	fscanf(promoter_f, "Contig=NT_%d", &g->NT);



	i = 0;
	c = ' ';
	while( int(c) != EOLN)
	{
		fscanf(promoter_f, "%c", &c);
		g->name[i++] = c;
		printf("GENENAME : %c    %d\n", g->name[i-1], int(g->name[i-1]));
	}
	g->name[i] = '\0';

	
	printf("Finished name fields Reading The gene ...\n");
	printf("Locus = %d Contig = %d gene name = %s\n",   g->locus, g->NT, g->name);
	
	c == ' ';




	i = 0;
	while( fscanf(promoter_f, "%c", &c) != EOF )
	{
		g->seq[i++] = c;
		
	}
	g->seq[i] = '\0';
*********************************************************************************************/


	fclose(expression_f);
	printf("Expressions were read successfully\n"); 

	return 0;

}






// Sort the genes such that the symbols and the names will match
long sort_all_genes(gene *g)
{
	long i, j;



//	char names[MAX_GENES_NUM][MAX_NAME_LENGTH];
//	char symbols[MAX_GENES_NUM][MAX_NAME_LENGTH];

	long indexes[MAX_GENES_NUM];


	long count = 0;


	char *temp_name;
	double *temp_exp;
	long temp_exp_len;
	int *temp_seq_exp;


	printf("TRUE %ld FALSE %ld\n", TRUE, FALSE);

	// start looking for identical <name, symbol> couples  
	for(i = 0; i < MAX_GENES_NUM-1; i++)
	{
		for(j = 0; j < EXP_GENES_NUM/*MAX_GENES_NUM*/; j++)
			if( strcmp(g[i].name, g[j].symbol ) == TRUE )
			{

				indexes[i] = j; 
			
				// Swap between them
				if( strcmp(g[j].name, g[j].symbol ) != TRUE ) // Not to 'Ruin' things which are the same..
				{
					count++;
					SWAP(g[i].symbol, g[j].symbol, temp_name);
					SWAP(g[i].exp_vec, g[j].exp_vec, temp_exp);
					SWAP(g[i].exp_len, g[j].exp_len, temp_exp_len);
				}
				break;


			}
		if( j == EXP_GENES_NUM)
			indexes[i] = -1;
	}

	printf("Found %ld Matches !!!\n", count);
	
	count = 0;
	for(i = 0; i < MAX_GENES_NUM-1; i++)
		if( strcmp(g[i].name, g[i].symbol ) == TRUE )
			count++;

	printf("FFFFFFF2Found %ld Matches !!!\n", count);






	// Now 'throw away' all that does not match ... 






	// Now we have to organize the genes according to what we have found
	for(i = 0; i < MAX_GENES_NUM-1; i++)
	{
		// check if we need to swap ...
		if( strcmp(g[i].name, g[i].symbol ) != TRUE )
		{

			j = i+1; 
			while( (strcmp(g[j].name, g[j].symbol ) != TRUE ) && (j < (MAX_GENES_NUM-1) ) )
				j++;

			if( /*strcmp(g[j].name, g[j].symbol ) == TRUE ) */ j < (MAX_GENES_NUM-1))
			{



				// Now do the Swap for all fields (maybe need a function for this ..)
				SWAP(g[i].symbol, g[j].symbol, temp_name);
				SWAP(g[i].exp_vec, g[j].exp_vec, temp_exp);
				SWAP(g[i].exp_len, g[j].exp_len, temp_exp_len);
				SWAP(g[i].accession, g[j].accession, temp_name);
				SWAP(g[i].afy_name, g[j].afy_name, temp_name);
				SWAP(g[i].locus, g[j].locus, temp_exp_len);
				SWAP(g[i].name, g[j].name, temp_name);
				SWAP(g[i].NT, g[j].NT, temp_exp_len);
				SWAP(g[i].seq, g[j].seq, temp_name);
				SWAP(g[i].seq_exp, g[j].seq_exp, temp_seq_exp);
				SWAP(g[i].seq_len, g[j].seq_len, temp_exp_len);


			}


		}
			
		

	}


	printf("Found %ld Matches !!!\n", count);
	count = 0;
	for(i = 0; i < MAX_GENES_NUM-1; i++)
		if( strcmp(g[i].name, g[i].symbol ) == TRUE )
			count++;
	printf("Now .. Found %ld Matches !!!\n", count);

	return count;


}



// What to do here ???? What does this function mean ? 
void search_for_motifs(gene *g)
{
}


// expand a sequence to contain all the overlapping.
// Are we sure that we want it ? what about the new function (for yuval ..)
long expand_seq(char *seq, long seq_len , int *seq_exp, long *seq_exp_len)
{
	long i, j;



	*seq_exp_len = seq_len - MAX_MODEL_SIZE+1;



	for(i = 0; i <= seq_len-MAX_MODEL_SIZE; i++)
	{
		seq_exp[i] = 0;

		for(j = 0; j < MAX_MODEL_SIZE; j++)
		{
			seq_exp[i] <<= 2;
			if(seq[i+j] == 'A')
				seq_exp[i] ^= 0;
			else if(seq[i+j] == 'T')
				seq_exp[i] ^= 1;
			else if(seq[i+j] == 'C')
				seq_exp[i] ^= 2;
			else if(seq[i+j] == 'G')
				seq_exp[i] ^= 3;
			else   // here not a legitamate letter. we 'through' it away
			{
				seq_exp[i] ^= 0x8000;
				break;
			}

		}



	}

	return 0;

}




// Calc correlations between genes and expressions ? 
long calc_correlation(gene *g, double *expression)
{
	return 0;
}


// Give the reverse compliment of a gene
void calc_reverse_complement(char *seq, long seq_len, char *comp_seq, long *comp_seq_len)
{
	long i;

	// Set length's  are equal
	*comp_seq_len = seq_len;

	// Calculate the Complement
	for(i = 0; i < seq_len; i++)
	{
		if(seq[seq_len - i-1] == 'A')
			comp_seq[i] = 'T';
		else if(seq[seq_len - i-1] == 'T')
			comp_seq[i] = 'A';
		else if(seq[seq_len - i-1] == 'C')
			comp_seq[i] = 'G';
		else if(seq[seq_len - i-1] == 'G')
			comp_seq[i] = 'C';
		else
			comp_seq[i] = seq[seq_len - i-1];
	}

	comp_seq[seq_len] = '\0';

}



// What the f*** is this ? 
long update_count(long *counts, char *seq, long seq_len)
{
	long i;
	long total = 0;
	long seq_exp_len;

	int *seq_exp = new int[seq_len];
	
	
	expand_seq(seq, seq_len , seq_exp, &seq_exp_len);



	for(i = 0; i < seq_exp_len; i++)
		if( (seq_exp[i]&0x8000) == 0)
		{
			counts[seq_exp[i]&0xFFF]++;
			total++;
		}
	
	delete seq_exp;

	return total;

}





// We currently assume a markov model of size up to 6
// This looks bad !!! The markov model table should  be always one (or two) dimensional, 
// regardless of the dimension of the markov model itself.
void markov_model_collect_statistics(char *seq, long seq_len,  long counts[4][4][4][4][4][4])
{
	long i;

	long i1, i2, i3, i4, i5, i6;


	// Make a 'naive' derichlet correction. Note - this is 'Big' 4^6 = 4096   
	for(i1 = 0; i1 < 4; i1++)
		for(i2 = 0; i2 < 4; i2++)
			for(i3 = 0; i3 < 4; i3++)
				for(i4 = 0; i4 < 4; i4++)
					for(i5 = 0; i5 < 4; i5++)
						for(i6 = 0; i6 < 4; i6++)
							counts[i1][i2][i3][i4][i5][i6] = 1;
	

	// Go over the sequence and collect the statistics
	for(i = 0; i < (seq_len - MAX_MODEL_SIZE + 1); i++)
		counts[seq[i]][seq[i+1]][seq[i+2]][seq[i+3]][seq[i+4]][seq[i+5]]++;

		
		
}




// Printing a probabilities array ..
void print_probs(double *probs, long dim)
{

	long i, j;



	printf("Probabilities :\n         A         T         C         G \n");
	printf("============================================\n");
	for(i = 0; i < (1UL << (2*dim-2)); i++)
	{
		for(j = 0; j < (dim-1); j++)
		{
			if( ((i >> (2*(dim-j-2)))&0x3) == 0 )
				printf("A");
			else if( ((i >> (2*(dim-j-2)))&0x3) == 1 )
				printf("T");
			else if( ((i >> (2*(dim-j-2)))&0x3) == 2 )
				printf("C");
			else if( ((i >> (2*(dim-j-2)))&0x3) == 3 )
				printf("G");
		
		}

		printf(" : ");
		for(j = 0; j < 4; j++)
			printf("%lf, ", probs[ (i << 2) ^ j]);
		printf("\n");
	}
	
	
	/**long i1, i2, i3, i4, i5, i6;
	
	for(i1 = 0; i1 < 4; i1++)
		for(i2 = 0; i2 < 4; i2++)
			for(i3 = 0; i3 < 4; i3++)
				for(i4 = 0; i4 < 4; i4++)
					for(i5 = 0; i5 < 4; i5++)
					{
						printf("[ %d %d %d %d %d ]  ", i1, i2, i3, i4, i5 );
						printf(" %lf, %lf, %lf, %lf\n", 
							probs[i1][i2][i3][i4][i5][0],
							probs[i1][i2][i3][i4][i5][1],
							probs[i1][i2][i3][i4][i5][2],
							probs[i1][i2][i3][i4][i5][3]);
					}

***/
}






// convert the counts into probabilities.
// derich is the amount of derichle correction
// dim is the dimention of the model required. (currently up to six ..)
void counts_to_probs(long *counts, long total_count, long dim, double derich, double *probs)
{
	long i, j;


	double total;



	// first "compress" the counts according to the dim
	for(i = 0; i < (1UL << (2*dim)); i++)
	{
		probs[i] = derich;
		for(j = 0; j < (1UL << (2*(MAX_MODEL_SIZE - dim))); j++)
		{
			probs[i] += counts[( i << (2*(MAX_MODEL_SIZE - dim)) ) + j];
		}


	}

	// calculate the total sum ...
	total = 0;
	for(i = 0; i < (1UL << (2*dim)); i++)
		total += probs[i];
	for(i = 0; i < (1UL << (2*dim)); i++)
		probs[i] /= total;


/***	
	// now turn counts into probabilities 
	for(i = 0; i < (1UL << (2*dim-2)); i++)
	{
		total = 0;
		for(j = 0; j < 4; j++)
			total += probs[ (i << 2) ^ j];
		for(j = 0; j < 4; j++)
			probs[ (i << 2) ^ j] /= total;

	}
***/

}


// Here we transfer from probabilities of cells to conditional probabilities.
void probs_to_markov_probs(double *probs, double *markov_probs, long dim)
{
	long i, j;

	double total;


	// first transfer the 'regular' probs into generation probs
	for(i = 0; i < (1UL << (2*dim-2)); i++)
	{
		total = 0;
		for(j = 0; j < 4; j++)
			total += probs[ (i << 2) ^ j];
		for(j = 0; j < 4; j++)
			markov_probs[ (i << 2) ^ j] = probs[ (i << 2) ^ j] / total;

	}
}


// Generate a random sequence according to the markov model 
// This is important for randomization !!!
void generate_random_seq(char *seq, long seq_len, double *probs, long dim)
{
	long i, j;

	double total;

	int cur_state, temp;

	double r;  // a random number between 0 to 1

	double gen_probs[1UL << 13]; // probabilities for generation of the model



	// first transfer the 'regular' probs into generation probs
	for(i = 0; i < (1UL << (2*dim-2)); i++)
	{
		total = 0;
		for(j = 0; j < 4; j++)
			total += probs[ (i << 2) ^ j];
		for(j = 0; j < 4; j++)
			gen_probs[ (i << 2) ^ j] = probs[ (i << 2) ^ j] / total;

	}

	cur_state = rand() &  ( (1UL << (2*dim)) - 1);


	// first do dome slak iterations to 'mix' everything, in order to get 
	// to the stationary distribution
	for(i = 0; i < 10/*SLAK_ITERS*/; i++)
	{
		temp = (cur_state << 2) & ( (1UL << (2*dim)) - 1);

		r = double(rand()) / double(RAND_MAX);

		if(r < gen_probs[temp])
			cur_state = temp;
		else if(r < gen_probs[temp] + gen_probs[temp^1])
			cur_state = temp^1;
		else if(r < gen_probs[temp] + gen_probs[temp^1] + gen_probs[temp^2])
			cur_state = temp^2;
		else
			cur_state = temp^3;



	}



	// Now generate the sequence ...
	for(i = 0; i < seq_len; i++)
	{

		temp = (cur_state << 2) & ( (1UL << (2*dim)) - 1);

		r = double(rand()) / double(RAND_MAX);



		if(r < gen_probs[temp])
		{
			cur_state = temp;
			seq[i] = 'A';
		}
		else if(r < gen_probs[temp] + gen_probs[temp^1])
		{
			cur_state = temp^1;
			seq[i] = 'T';
		}
		else if(r < gen_probs[temp] + gen_probs[temp^1] + gen_probs[temp^2])
		{
			cur_state = temp^2;
			seq[i] = 'C';
		}
		else
		{
			cur_state = temp^3;
			seq[i] = 'G';
		}

	}



}




// Generate a random sequence according to the markov model 
// This is important for randomization !!!
void generate_random_seq_word(word *seq, long seq_len, double *probs, long dim)
{
	long i, j;

	double total;

	int cur_state, temp;

	double r;  // a random number between 0 to 1

	double gen_probs[1UL << 13]; // probabilities for generation of the model



	// first transfer the 'regular' probs into generation probs
	for(i = 0; i < (1UL << (2*dim-2)); i++)
	{
		total = 0;
		for(j = 0; j < 4; j++)
			total += probs[ (i << 2) ^ j];
		for(j = 0; j < 4; j++)
			gen_probs[ (i << 2) ^ j] = probs[ (i << 2) ^ j] / total;

	}

	cur_state = rand() &  ( (1UL << (2*dim)) - 1);


	// first do dome slak iterations to 'mix' everything, in order to get 
	// to the stationary distribution
	for(i = 0; i < 10/*SLAK_ITERS*/; i++)
	{
		temp = (cur_state << 2) & ( (1UL << (2*dim)) - 1);

		r = double(rand()) / double(RAND_MAX);

		if(r < gen_probs[temp])
			cur_state = temp;
		else if(r < gen_probs[temp] + gen_probs[temp^1])
			cur_state = temp^1;
		else if(r < gen_probs[temp] + gen_probs[temp^1] + gen_probs[temp^2])
			cur_state = temp^2;
		else
			cur_state = temp^3;



	}


	// Now generate the sequence ...
	for(i = 0; i < seq_len; i++)
	{

		temp = (cur_state << 2) & ( (1UL << (2*dim)) - 1);

		r = double(rand()) / double(RAND_MAX);

		if(r < gen_probs[temp])
		{
			cur_state = temp;
			seq[i/HALF_WORD_SIZE] ^= 0; // 'A'
		}
		else if(r < gen_probs[temp] + gen_probs[temp^1])
		{
			cur_state = temp^1;
			seq[i/HALF_WORD_SIZE] ^= (0x1 << (2* (i&0x7)));  // 'T' We start from least to most

		}
		else if(r < gen_probs[temp] + gen_probs[temp^1] + gen_probs[temp^2])
		{
			cur_state = temp^2;
			seq[i/HALF_WORD_SIZE] ^= (0x2 << (2* (i&0x7)));  // 'C' 
		}
		else
		{
			cur_state = temp^3;
			seq[i/HALF_WORD_SIZE] ^= (0x3 << (2* (i&0x7)));  // 'G' 
		}

	}


}



// This transfers the probabilities of the counts into weights used to give a score, 
// by applying log transform
void probs_to_weights(double *probs, long dim, double *weights)
{

	long i, j;

	double total;

	// transfer the 'regular' probs into generation probs, and then log transform them
	for(i = 0; i < (1UL << (2*dim-2)); i++)
	{
		total = 0;
		for(j = 0; j < 4; j++)
			total += probs[ (i << 2) ^ j];
		for(j = 0; j < 4; j++)
			weights[ (i << 2) ^ j] = log( probs[ (i << 2) ^ j] / total );
	}


}



// Counts of the markov model to probabilities - in the old style ...
void count_to_probs(long counts[4][4][4][4][4][4], double probs[4][4][4][4][4][4])
{
	long i1, i2, i3, i4, i5;
	
	double total;

	for(i1 = 0; i1 < 4; i1++)
		for(i2 = 0; i2 < 4; i2++)
			for(i3 = 0; i3 < 4; i3++)
				for(i4 = 0; i4 < 4; i4++)
					for(i5 = 0; i5 < 4; i5++)
					{
						total = double( counts[i1][i2][i3][i4][i5][0] + 
										counts[i1][i2][i3][i4][i5][1] + 
										counts[i1][i2][i3][i4][i5][2] + 
										counts[i1][i2][i3][i4][i5][3] );

						probs[i1][i2][i3][i4][i5][0] = counts[i1][i2][i3][i4][i5][0] / total;
						probs[i1][i2][i3][i4][i5][1] = counts[i1][i2][i3][i4][i5][1] / total;
						probs[i1][i2][i3][i4][i5][2] = counts[i1][i2][i3][i4][i5][2] / total;
						probs[i1][i2][i3][i4][i5][3] = counts[i1][i2][i3][i4][i5][3] / total;


					}
	
	

}


// Read matrices from files
void read_TF_matrices(FILE *tfname_f, char *tfname_p, FILE *gename_f, char *gename_p, 
					  long pssm[L_TF][4][MAX_L], long pssm_len[MAX_L], 
					  char TF_NAMES[L_TF][MAX_NAME_LENGTH], char gene_names[L_GENES][MAX_NAME_LENGTH])
{

	long i, j, k;

	char c;





/***/
	// open the file to enable reading from it the gene names
	gename_f = fopen(gename_p, "r");

	for(i = 0; i < L_GENES; i++)
		fscanf(gename_f, "%s\n", gene_names[i]);

	// close file
	fclose(gename_f);

/***/

	// open the file to enable reading from it the TF names
	tfname_f = fopen(tfname_p, "r");


	for(i = 0; i < L_TF; i++)
		fscanf(tfname_f, "%s\n", TF_NAMES[i]);

	// close file
	fclose(tfname_f);


	// finished reading names .....


	FILE *pssm_f;          // pssm of TF file 
//	char *pssm_p = "SC_gene_results.txt";   


	for(i = 0; i < L_TF; i++)
	{
		// open the file to enable reading from it
		pssm_f = fopen(TF_NAMES[i], "r");

		j = 0;

		c = 'f';
		while( int(c) != EOLN)
		{
			fscanf(pssm_f, "%ld", &(pssm[i][0][j++]));
			fscanf(pssm_f, "%c", &c);
		}


		pssm_len[i] = j;
		printf("TF NAME : %s PSSM LEN : %ld\n", TF_NAMES[i], j);

		for(k = 1; k < 4; k++)
			for(j = 0; j < pssm_len[i]; j++)
				fscanf(pssm_f, "%ld", &(pssm[i][k][j]));
		


		fclose(pssm_f);

	}






	

}



// Read the lengthes of the genes
void read_lengths(FILE *lengths_f, char *lengths_p, long genes_lengths[L_GENES])
{
	long i;

	// open the file to enable reading from it
	lengths_f = fopen(lengths_p, "r");

	for(i = 0; i < L_GENES; i++)
		fscanf(lengths_f, "%ld", &(genes_lengths[i]));

	// close file
	fclose(lengths_f);


}



// Read thresholds to be used
void read_threshes(FILE *threshs_f, char *threshs_p, double genes_threshes[L_TF][L_GENES])
{
	long i, j;


	// open the file to enable reading from it
	threshs_f = fopen(threshs_p, "r");

	for(i = 0; i < L_TF; i++)
		for(j = 0; j < L_GENES; j++)
			fscanf(threshs_f, "%lf", &(genes_threshes[i][j]));

	// close file
	fclose(threshs_f);

}


// Read markov probs from a file 
void read_markov_probs(FILE *probs_f, char *probs_p, double *probs)
{
	long i;		
	
	// open the file to enable reading from it
	probs_f = fopen(probs_p, "r");

	// now read the probabilities
	for(i = 0; i < 16; i++)
		fscanf(probs_f, "%lf", &(probs[i]));	

	// close file
	fclose(probs_f);
}




// Write markov probs from a file 
void write_probs(FILE *probs_f, char *probs_p, double probs[L_TF][L_GENES])
{
	long i, j;
			
	// open the file to enable reading from it
	probs_f = fopen(probs_p, "w");

	// now read the probabilities

	for(i = 0; i < L_GENES; i++)
	{
		for(j = 0; j < L_TF; j++)
			fprintf(probs_f, "%lf ", probs[j][i]);	
		fprintf(probs_f, "\n");
	}


	// close file
	fclose(probs_f);
}




// Transfer the position specific score matrices and the background model 
// into weights which are their logs, in order to use the weights for easy scoring 
// of potential binding sites.
void pssm_and_background_to_weights(long pssm[4][MAX_L], long pssm_len, double derich, double *back_probs, 
									double w[4][4][MAX_L])
{
	long i, j, k;
	
	double pssm_w[4][MAX_L];	
	double total;

	
	// first generate the pssm weights
	for(i = 0; i < pssm_len; i++)
	{
		total = pssm[0][i] + pssm[1][i] + pssm[2][i] + pssm[3][i] + 4*derich;

		for(j = 0; j < 4; j++)
			pssm_w[j][i] = log(  (pssm[j][i] + derich)/total );

	}


	// Now take all into account
	for(i = 0; i < pssm_len; i++)
		for(j = 0; j < 4; j++)
			for(k = 0; k < 4; k++)
				w[j][k][i] = pssm_w[k/*j*/][i] - /*1*/1*log(back_probs[/*4*k+j*/4*j+k]);

}



// calculate the score of a given sequence according to a site 
double calc_site_score(double w[4][4][MAX_L], long l, long seq[MAX_L])
{
	long i;

	double score = 0;

	for(i = 0; i < l; i++)
		score += w[ seq[i] ] [ seq[i+1] ] [i];

	return score;

}



///////////////////////////////////////////////////////////////////////////////////////
/// New Functions For Yuval's Project - 11.11.2003
///////////////////////////////////////////////////////////////////////////////////////


/****** put it in a  mex file 
// len is in nucleotides 
void print_packed_dna(word *seqs, long len)
{
	long i, j;
	word tmp;

	for(i = 0; i < len; i++)
	{

		tmp = seqs[i/16];
		tmp = tmp >> (2*(i&0xF));
		tmp = tmp&0x3;
		tmp++;

//		tmp = ((seqs[i/16] << (2* (i&0xF)))&0x3) + 1;


//		tmp = ((seqs[i/16] >> ((i&0xF)*2))&0x3) + 1;
		printf("%lx ", tmp);
	}
	printf("\n");

}
*******/

// Calculate the scores corresponding to position weight matrices, for many sites 
// We assume that seqs is one big array ... (need to take care of that in matlab ..)
// We assume that the scores array is already allocated. 
long calc_all_sites_scores(double weights[4][MAX_L], long L, word *seqs[MAX_NUM_SEQS], 
						   double *scores[MAX_NUM_SEQS], long num_seqs, long *seqs_lens)
{
	long i, j, g;	

	// Change to avoid checking the edges (Need to change back at the end !!!)
	for(g = 0; g < num_seqs; g++)
		seqs_lens[g] = seqs_lens[g] - L + 1;

	// Initilize all the scores to zero
	for(g = 0; g < num_seqs; g++)
		for(j = 0; j < seqs_lens[g]; j++)
			scores[g][j] = 0;

	// score the binding sites	
	for(g = 0; g < num_seqs; g++)
		for(j = 0; j < L; j++)
			for(i = 0; i < seqs_lens[g]; i++)
				scores[g][i] += weights[(seqs[g][(i+j)/HALF_WORD_SIZE] >> (2*   ((i+j)&0xF)  ))&0x3][j];

	// Need to change back at the end !!!
	for(g = 0; g < num_seqs; g++)
		seqs_lens[g] = seqs_lens[g] + L - 1;
	return 0;
}


// Calculate the scores corresponding to position weight matrices, for many sites 
// We assume that seqs is one big array ... (need to take care of that in matlab ..)
// We assume that the scores array is already allocated. 
// Essentialy the same function as calc_all_sites_scores but for singles
long calc_all_sites_scores_single(float weights[4][MAX_L], long L, word *seqs[MAX_NUM_SEQS], 
						   float *scores[MAX_NUM_SEQS], long num_seqs, long *seqs_lens)
{
	long i, j, g;	

	// Change to avoid checking the edges (Need to change back at the end !!!)
	for(g = 0; g < num_seqs; g++)
		seqs_lens[g] = seqs_lens[g] - L + 1;

	// Initilize all the scores to zero
	for(g = 0; g < num_seqs; g++)
		for(j = 0; j < seqs_lens[g]; j++)
			scores[g][j] = 0;

	// score the binding sites	(assumption: L <= 16)
	for(g = 0; g < num_seqs; g++)
		for(j = 0; j < L; j++)
			for(i = 0; i < seqs_lens[g]; i++)
				scores[g][i] += weights[(seqs[g][(i+j)/HALF_WORD_SIZE] >> (2*   ((i+j)&0xF)  ))&0x3][j];

	// Need to change back at the end !!!
	for(g = 0; g < num_seqs; g++)
		seqs_lens[g] = seqs_lens[g] + L - 1;
	return 0;
}







// Get an array of scores, and output the locations & scores of the ones which have past the threshold.
// In the future it is possible to let the routine select it's own threshold in some smart way.
long get_scores_above_threshold(long is_double, double *scores[MAX_NUM_SEQS], float *scores_single[MAX_NUM_SEQS], 
								long num_seqs, long *seqs_lens, double threshold, float threshold_single,
								long remove_overlap_distance, long smart_thresh, 
								long *out_x, long *out_y, 
								double *out_scores, float *out_scores_single) // Those are allocated OUTSIDE the function
{
	long i, j, k, g;	
	long dist = remove_overlap_distance; 
	long count_above = 0;
	long counter = 0;
	double local_threshold = threshold;
	float local_threshold_single = threshold_single; 

//	long *indexes; 

	// Now determine the threshold if needed. Here at the beginning threshold is the FRACTION !!! 
//		if(smart_threshold == 1)
//		for(g = 0; g < num_seqs; g++)
//		{
//			DoQuicksort(scores[g], seqs_lens[g], indexes)
//		}
			
	if(is_double)
	{
		for(g = 0; g < num_seqs; g++) 		// First go over all scores to see how many are above threshold
			for(j = 0; j < seqs_lens[g]; j++)
				if(scores[g][j] >= local_threshold)
					count_above++;
		for(g = 0; g < num_seqs; g++) 			// Now loop again over all scores to fill the new arrays
			for(j = 0; j < seqs_lens[g]; j++)
				if(scores[g][j] >= local_threshold)
				{
					out_x[counter] = g+1; 
					out_y[counter] = j+1;   // we do plus 1 because we go to matlab ...
					out_scores[counter++] = scores[g][j];
				}
	}
	else
	{		
		for(g = 0; g < num_seqs; g++)
			for(j = 0; j < seqs_lens[g]; j++)
				if(scores_single[g][j] >= local_threshold_single)
					count_above++;
		for(g = 0; g < num_seqs; g++) 			// Now loop again over all scores to fill the new arrays
			for(j = 0; j < seqs_lens[g]; j++)
				if(scores_single[g][j] >= local_threshold_single)
				{
					out_x[counter] = g+1; 
					out_y[counter] = j+1;   // we do plus 1 because we go to matlab ...
					out_scores_single[counter++] = scores_single[g][j];
				}
	}
			
	long count_good = count_above; // The number of good scores left
	if(dist > (-1)) // -1 means no overlap removal
	{
		long *good_ind = new long[count_above];
		for(i = 0; i < count_above; i++)
			good_ind[i] = 1;
		for(i = 0; i < count_above; i++)
		{
			for(k = i+1; k < count_above; k++)  // go forward
			{
				if(out_x[k] != out_x[i])  // we got to some other gene
					break; 
				if(out_y[k] > out_y[i] + dist)  // here we are more than dist away
					break;
				if(is_double)
				{
					if(out_scores[k] > out_scores[i])  // here we got a better scores !!!
					{
						good_ind[i] = 0; count_good--;
						break;
					}
				}
				else
				{
					if(out_scores_single[k] > out_scores_single[i])  // here we got a better scores !!!
					{
						good_ind[i] = 0; count_good--;
						break;
					}
				}
			}
			if(good_ind[i])
				for(k = i-1; k >=0; k--)  // go also backwards if needed
				{
					if(out_x[k] != out_x[i])  // we got to some other gene
						break; 
					if(out_y[k] < out_y[i] - dist)  // here we are more than dist away
						break;
					if(is_double)
					{
						if(out_scores[k] > out_scores[i])  // here we got a better scores !!!
						{
							good_ind[i] = 0; count_good--;
							break;
						}
					}
					else
					{
						if(out_scores_single[k] > out_scores_single[i])  // here we got a better scores !!!
						{
							good_ind[i] = 0; count_good--;
							break;
						}
					}
				}
		}

		// Now that we know the good indices, we can set the array.
		long *new_x = new long[count_good];
		long *new_y = new long[count_good];
		double *new_scores  = new double[count_good];
		float *new_scores_single  = new float[count_good];
		/**
		if(is_double)
			new_scores = new double[count_good];
		else
			new_scores_single = new float[count_good];
		**/

		counter = 0;
		for(i = 0; i < count_above; i++)
			if(good_ind[i])
			{
				new_x[counter] = out_x[i];
				new_y[counter] = out_y[i];
				if(is_double)
					new_scores[counter++] = out_scores[i];
				else
					new_scores_single[counter++] = out_scores_single[i];
			}

		// Set the pointers back  - We cannot set pointers since we will 'lose' them, so we copy whole the values
		for(i = 0; i < count_good; i++)
		{
			out_x[i] = new_x[i];
			out_y[i] = new_y[i];
			if(is_double)
				out_scores[i] = new_scores[i];
			else
				out_scores_single[i] = new_scores_single[i];
		}
		delete good_ind;
		delete new_x;
		delete new_y;
		delete new_scores;
		delete new_scores_single;
	}
	return count_good;
}





// Spread the sequences in order to make the binding sites score calculations easier
// We assume that spreaded_seqs are already allocated with enough memory. 
// The preciese memory needed is : seqs_len * num_genes * 2 bits ...
// THIS IS NOT NEEDED RIGHT NOW !!!
long spread_sequences_to_binding_sites(word *seqs[MAX_SEQ_LEN], long num_genes, long seqs_len, word *spreaded_seqs, long L)
{

	long i, j;


	// read and spread gene by gene or location by location ...
	for(i = 0; i < num_genes; i++)
	{
		// spread the following gene
		for(j = 0; j < seqs_len - L+1; j++)
		{
			// The letters are aligned to the least significant bit ! 
			spreaded_seqs[j] = (seqs[i][j/HALF_WORD_SIZE] >> (2*(j&0xF))) ^ 
							   (seqs[i][(j+1)/HALF_WORD_SIZE] >> 9); // ??? This should change !!
		}

	}

	return 0;

}


///////////////////////////////////////////////////////////////////////////////////////
/// End New Functions For Yuval's Project - 11.11.2003
///////////////////////////////////////////////////////////////////////////////////////




int main()
{
	long i, j, l, L;	
		
	char my_seq[SEQ_LEN+1];
	char my_comp_seq[SEQ_LEN+1];
	char *gen_seq;



	long my_seq_len, my_comp_seq_len, gen_seq_len;

	long my_seq_exp_len;

	long temp;

	long my_dim = 2;
	double my_derich = 0.5;


	long my_counts[1UL << 13];
	long gen_counts[1UL << 13];
	
	
	double my_probs[1UL << 13];
	double my_copy_probs[1UL << 13];
	double gen_probs[1UL << 13];
//	double markov_probs[1UL << 13];
	
	
	double *my_dist;

	long my_total;
	long gen_total;

//	long *exp_mat[MAX_SAMPLES];

//	exp_vec gev[MAX_GENES_NUM];	


	gene g;

	gene genes[MAX_GENES_NUM];

//	gene gev[7000];
//	char labels[MAX_SAMPLES][MAX_NAME_LENGTH];


	FILE *promoter_f;		// promoter file of the gene 
	FILE *human13k_f;		// promoter file of all 13k human genes name
//	FILE *expression_f;		// expression file of all genes name
	FILE *probs_f;          // probabilities calculated
	FILE *markov_probs_f;          // probabilities calculated
	FILE *threshs_f;          // thresholds calculated
	FILE *lengths_f;          // lengths of genes 
	FILE *tfname_f;          // TF names
	FILE *gename_f;          // gene names
	FILE *pvalues_f;          // p values for TF's and genes
	FILE *libipvalues_f;      // p values from libi (assumes only singeltons statistics)

	char *promoter_p = "first_gene.txt";			// 	promoter file name
	char *human13k_p = "HumanPromoters13K.txt"; //"few_genes.txt";		// promoter file of all 13k human genes name
	char *expression_p = "new_exp.txt";					// expression file of all genes name
	char *probs_p = "probs.txt";
	char *markov_probs_p = "SC_gene_names_markov_1.txt";  

	char *lengths_p = "SC_gene_names_lens.txt";							// lengths of genes 
	char *tfname_p = "SC_TF_names.txt";          // TF names
	char *gename_p = "SC_gene_names.txt";          // gene names

	char *libipvalues_p = "get_TFs_p_val_script_SC_gene_names_2_non_uni_1_or_1_4.txt"; 
	//"SC_TF_genes_p_values.txt";      // p values from libi (assumes only singeltons statistics)

#ifdef DEBUG
	char *threshs_p = "SC_genes_scores_indep_1_4.txt";   // scores (log prob) calculated
	char *pvalues_p = "SC_TF_genes_p_values_DEBUG.txt";   // p values for TF's and genes
#else
	char *threshs_p = "SC_gene_results.txt";   // scores (log prob) calculated
	char *pvalues_p = "SC_TF_genes_p_values.txt";   // p values for TF's and genes
#endif // DEBUG

	printf("start ,,,\n");

	// get the gene's length
	g.seq_len = find_gene_length(promoter_f, promoter_p);


	// allocating memory for the gene
	g.seq = new char[g.seq_len];
	g.name = new char[MAX_NAME_LENGTH];


	// first read gene to deal with 
	read_gene(promoter_f, promoter_p, &g);



	//srand(sysclock);

	// Init the sequence 
	my_seq_len = SEQ_LEN;


	// Set random DNA
	for(i = 0; i < my_seq_len; i++)
	{	
		temp = rand();
			

		if( (temp&0x3) == 0)
			my_seq[i] = 'A';
		else 	if( (temp&0x3) == 1)
			my_seq[i] = 'T';
		else 	if( (temp&0x3) == 2)
			my_seq[i] = 'C';
		else 	
			my_seq[i] = 'G';

	//		my_seq[i] = temp&0x3;
	}
	my_seq[my_seq_len] = '\0';


	calc_reverse_complement(my_seq, my_seq_len, my_comp_seq, &my_comp_seq_len);

	print_dna(my_seq, 10);
	print_dna(my_comp_seq, 10);


	g.seq_exp = new int[g.seq_len];
	expand_seq(g.seq, g.seq_len , g.seq_exp, &my_seq_exp_len);
		




//////////////////////////////////////////////////////////////////////////////////////


	// Now work on all genes ...
	for(i = 0; i < MAX_GENES_NUM; i++)
	{
		genes[i].seq = new char[1220];
		genes[i].name = new char[MAX_NAME_LENGTH];
		genes[i].symbol = new char[MAX_NAME_LENGTH];
		genes[i].afy_name = new char[MAX_NAME_LENGTH];
		genes[i].exp_vec = new double[RAMAS_SAMPLES];
		genes[i].seq_len = 1201;
		genes[i].exp_len = RAMAS_SAMPLES-2;
	}





	// first read gene to deal with 
	read_all_genes(human13k_f, human13k_p, genes);



	// Now update our statistics ....
	for(i = 0; i < (1UL << 13); i++)
		my_counts[i] = 0;

	// do counts
	for(i = 0; i < MAX_GENES_NUM; i++)
		my_total += update_count(my_counts, genes[i].seq, genes[i].seq_len);

	counts_to_probs(my_counts, my_total, my_dim,  my_derich, my_probs);

	printf("collected probs :\n");
	print_probs(my_probs, my_dim);


	// open the file to enable writing to it
	probs_f = fopen(probs_p, "w+");





#ifdef DO_SORT
	read_expression(expression_f, expression_p, genes, labels );

	// print selected genes
	sort_all_genes(genes);
	printf("sorted genes ..\n");

#endif




	// New 11.12.2003 - We do not need anything below ...
//	exit(99999); 
//////////////////////////////////////////////////////////////////////////////////////////////////


/// Here do the part of the TF's stuff

	long STAT[MAX_L][4] = {
		{ 7,      5,      9,      7, },     
		{ 7,      7,      8,      6, },    
		{ 12,      6,      5,     7, },   
		{ 5,      7,      9,      9, },    
		{ 9,      5,      6,     10, },     
		{ 2,      0,      0,     28, },    
		{ 0,      0,      0,     30, },     
		{ 0,      0,      0,     30, },     
		{ 0,     30,      0,      0, },     
		{ 2,     21,      4,      3, },     
		{ 4,     14,     11,      1, },     
		{ 3,      4,     20,      3, },     
		{ 0,      0,     30,      0, },     
		{ 30,      0,      0,     0, },     
		{ 30,      0,      0,     0, },     
		{ 30,      0,      0,     0, },     
		{ 4,      8,      5,     12, },     
		{ 7,      6,      9,      6, },     
		{ 5,      8,      8,      5, },    
		{ 6,      8,      2,     10, },    
		{ 9,      3,      7,      7, }};


//	double site_probs[4][4][MAX_L];
	double w[4][4][MAX_L]; 	
	long m[4][4][MAX_L];
	double my_res = 0.001;
	long my_thresh = 200;


	l = 21;

	// Do derichle correction
	for(i = 0; i < 4; i++)
		for(j = 0; j < 21; j++)
#define DO_EXAMPLE
#ifdef DO_EXAMPLE
			STAT[j][i]++;

//	probs_to_markov_probs(my_probs, markov_probs, my_dim);
#endif

	return 0; 

}





































