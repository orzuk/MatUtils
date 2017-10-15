// File general.h - Macros, prototypes and constant of the Rijndael algorithm


#ifndef _GENERAL_H_
#define _GENERAL_H_


// Include files needed
#include <string.h>
#include <stdio.h>
#include <iostream.h>
#include <stdlib.h>
#include <time.h>




// New types defined 
typedef unsigned long int word;   // machine word (32 bit on windows - pc)

typedef unsigned char byte;       // a byte

// Macroes 
#define BIT(block, num)  ( ((block) >> (num)) & (1UL) )  // get the i-th bit from a word


// inverse the order of the bits in a word
#define INVERSE(block) \
	(block) = (((block)&0x0000FFFF) << 16 ) ^ (((block) >> 16)&0x0000FFFF); \
	(block) = (((block)&0x00FF00FF) << 8 ) ^ (((block) >> 8)&0x00FF00FF); \
	(block) = (((block)&0x0F0F0F0F) << 4 ) ^ (((block) >> 4)&0x0F0F0F0F); \
	(block) = (((block)&0x33333333) << 2 ) ^ (((block) >> 2)&0x33333333); \
	(block) = (((block)&0x55555555) << 1 ) ^ (((block) >> 1)&0x55555555); 

// inverse the even and the odd bits in a block
#define INV_PAR(block) \
	(block) = (((block)&0x55555555) << 1 ) ^ (((block) >> 1)&0x55555555); 

// inverse the order of bytes in a word (this is for changing indians)
#define INDY(block) \
	(block) = (((block)&0x0000FFFF) << 16 ) ^ (((block) >> 16)&0x0000FFFF); \
	(block) = (((block)&0x00FF00FF) << 8 ) ^ (((block) >> 8)&0x00FF00FF); 


#define XTIME(byte) ((byte) << 1) ^ BIT((byte), 7) * (0x1B)   // Multiply a byte by x 


#define SWAP_S(addr) \
	(addr) = ((addr)&0x20) ^ (((addr)&0x1) << 4) ^ (((addr)&0x1E) >> 1);

#define MAX(x,y) (x) > (y) ? (x) : (y)
#define MIN(x,y) (x) < (y) ? (x) : (y)


// Global constants 
#define LENGTH0 4   // The length of rijndael key and block. (versions with other key 
#define LENGTH1 6					  // lengthes are not supported yet)
#define LENGTH2 8

#define WORDLEN 32 // length of Rijndael (and machine) word

#define MAX_ROUNDS (LENGTH2 + 6)   // Maximum possible Rijndael rounds (14)

#define MAX_SHIFT 4				   // Maximum possible Rijndael shift

#define ENCRYPT 0
#define DECRYPT 1

#define TRUE 1 
#define FALSE 0


// machine modes 
#define ECB 0   // Electronic Code Book
#define CBC 1	// Cipher Block Chaining
#define OFB 2	// Output Feedback
#define CFB 3	// Cipher Feedback




// levels of printing
#define NO_PRINT    (-1)
#define PRINT_LEVEL_0 0
#define PRINT_LEVEL_1 1
#define PRINT_LEVEL_2 2
#define PRINT_LEVEL_3 3

// The indians of the machine 
#define LITTLE  0
#define BIG     1


//#define INDIANS  BIG 
#define INDIANS  LITTLE 


// The sbox of the Byte Sub
static byte sbox[256] = { 
 99, 124, 119, 123, 242, 107, 111, 197,  48,   1, 103,  43, 254, 215, 171, 118, 
202, 130, 201, 125, 250,  89,  71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 
183, 253, 147,  38,  54,  63, 247, 204,  52, 165, 229, 241, 113, 216,  49,  21, 
  4, 199,  35, 195,  24, 150,   5, 154,   7,  18, 128, 226, 235,  39, 178, 117, 
  9, 131,  44,  26,  27, 110,  90, 160,  82,  59, 214, 179,  41, 227,  47, 132, 
 83, 209,   0, 237,  32, 252, 177,  91, 106, 203, 190,  57,  74,  76,  88, 207, 
208, 239, 170, 251,  67,  77,  51, 133,  69, 249,   2, 127,  80,  60, 159, 168, 
 81, 163,  64, 143, 146, 157,  56, 245, 188, 182, 218,  33,  16, 255, 243, 210, 
205,  12,  19, 236,  95, 151,  68,  23, 196, 167, 126,  61, 100,  93,  25, 115, 
 96, 129,  79, 220,  34,  42, 144, 136,  70, 238, 184,  20, 222,  94,  11, 219, 
224,  50,  58,  10,  73,   6,  36,  92, 194, 211, 172,  98, 145, 149, 228, 121, 
231, 200,  55, 109, 141, 213,  78, 169, 108,  86, 244, 234, 101, 122, 174,   8, 
186, 120,  37,  46,  28, 166, 180, 198, 232, 221, 116,  31,  75, 189, 139, 138, 
112,  62, 181, 102,  72,   3, 246,  14,  97,  53,  87, 185, 134, 193,  29, 158, 
225, 248, 152,  17, 105, 217, 142, 148, 155,  30, 135, 233, 206,  85,  40, 223, 
140, 161, 137,  13, 191, 230,  66, 104,  65, 153,  45,  15, 176,  84, 187,  22, 
};

// The inverse sbox of the inverse Byte Sub
static byte inv_sbox[256] = {
 82,   9, 106, 213,  48,  54, 165,  56, 191,  64, 163, 158, 129, 243, 215, 251, 
124, 227,  57, 130, 155,  47, 255, 135,  52, 142,  67,  68, 196, 222, 233, 203, 
 84, 123, 148,  50, 166, 194,  35,  61, 238,  76, 149,  11,  66, 250, 195,  78, 
  8,  46, 161, 102,  40, 217,  36, 178, 118,  91, 162,  73, 109, 139, 209,  37, 
114, 248, 246, 100, 134, 104, 152,  22, 212, 164,  92, 204,  93, 101, 182, 146, 
108, 112,  72,  80, 253, 237, 185, 218,  94,  21,  70,  87, 167, 141, 157, 132, 
144, 216, 171,   0, 140, 188, 211,  10, 247, 228,  88,   5, 184, 179,  69,   6, 
208,  44,  30, 143, 202,  63,  15,   2, 193, 175, 189,   3,   1,  19, 138, 107, 
 58, 145,  17,  65,  79, 103, 220, 234, 151, 242, 207, 206, 240, 180, 230, 115, 
150, 172, 116,  34, 231, 173,  53, 133, 226, 249,  55, 232,  28, 117, 223, 110, 
 71, 241,  26, 113,  29,  41, 197, 137, 111, 183,  98,  14, 170,  24, 190,  27, 
252,  86,  62,  75, 198, 210, 121,  32, 154, 219, 192, 254, 120, 205,  90, 244, 
 31, 221, 168,  51, 136,   7, 199,  49, 177,  18,  16,  89,  39, 128, 236,  95, 
 96,  81, 127, 169,  25, 181,  74,  13,  45, 229, 122, 159, 147, 201, 156, 239, 
160, 224,  59,  77, 174,  42, 245, 176, 200, 235, 187,  60, 131,  83, 153,  97, 
 23,  43,   4, 126, 186, 119, 214,  38, 225, 105,  20,  99,  85,  33,  12, 125, 
};




static long pop_tab[16] = {
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4
};     // the pop count of a nibble


// The round constant used for the key expantion
static byte round_const[(MAX_ROUNDS+1)*2] = { 
	0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a, 
	0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91, };




// The machine struct. This contains all the neccessary information needed
// on the Rijndael machine.
typedef struct {


long blocklength;   // length of rijndael block
long keylength;     // length of rijndael key
long rounds;        // number of rijndael rounds to preform  
long round;       // what round of Rijndael is performed now 
long mode;        // 0 - ECB 1 - CBC 2 - OFB 3 - CFB
long flag;        // 0 - encrypt 1 - decrypt
long use_data;    // 1 - use the data block in the parms file, 0 - use only plain file
long pflag;       // flag saying the level of printing
long shifts[4];   // by how much to shifts the rows in the 
				  // shiftrow transformation
long expanded;    // flag saying if the key expantion was already done or not

byte round_const[(MAX_ROUNDS+1)*2];  // The round constants for the key expantion

byte sbox[256];   // The sbox of the bytesub 
byte inv_sbox[256];  // The inverse of the sbox

byte state[LENGTH2][4];  // the state array of bytes
byte key[LENGTH2][4];    // the key array of bytes
byte round_key[MAX_ROUNDS+1][LENGTH2][4];   // the keys that are used in each round
byte inv_round_key[MAX_ROUNDS+1][LENGTH2][4];   // the keys that are used in 
												// each round in the inverse cipher 


word *w_state;
word *w_key;
word *w_round_key[MAX_ROUNDS+1];	// These are word pointers that refer to the 
									// same data that is in the bytes.

word iv[LENGTH2][4];         // initialization vector , needed for most modes.
word initial_iv[LENGTH2][4];   // the initial initialization vector



long datalength;  // the length of data to encrypt (in 32 bit blocks)



} mach_state;





// Functions prototypes

long mix_colomn(mach_state *mt);
long shift_row(mach_state *mt);
long byte_sub(mach_state *mt);
long add_round_key(mach_state *mt);
long perform_round(mach_state *mt);
long perform_rijndael(mach_state *mt);
long expand_key(mach_state *mt);
long initilize_machine(mach_state * mt);
long print_state(mach_state *mt);
long master_rijndael(word *data_in, mach_state *mt, FILE *pr_f,
						word *data_out);
long load_state(mach_state *mt, word *data_in);
long save_state(mach_state *mt, word *data_out);


// The inverse Rijndael functions 
long inv_mix_colomn(mach_state *mt);
long inv_shift_row(mach_state *mt);
long inv_byte_sub(mach_state *mt);
long add_inv_round_key(mach_state *mt);
long perform_inv_round(mach_state *mt);
long perform_inv_rijndael(mach_state *mt);
long inv_mix_column_key(mach_state *mt);
long inv_expand_key(mach_state *mt);

// General functions
long read_parameters(FILE *parm_f, char *parm_p, mach_state * mt);
long read_data(FILE *data_in_f, char *data_in_p, long datalength,
			   word *data_in);
long write_data(FILE *data_out_f, char *data_out_p, long datalength,
			   word *data_out);
long calc_data_length(FILE *data_in_f, char *data_in_p);

long x_inv_box(mach_state *mt);
byte poppar(byte block);

#endif // GENERAL_H 