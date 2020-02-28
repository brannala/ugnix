#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <gmodule.h>
#include <math.h>
#include <unistd.h>
#include <stdbool.h>
#include <regex.h>

#define MAXLN 100
#define MAXALLNM 50
#define MAXALL 500
#define MAXIND 10000
#define MAXPOP 500
#define MAXLOCI 100000

/* constants */

#define PLL_STRING(x) #x
#define PLL_C2S(x) PLL_STRING(x)

#define VERSION_MAJOR 0
#define VERSION_MINOR 2
#define VERSION_PATCH 1

#define PROG_VERSION "v" PLL_C2S(VERSION_MAJOR) "." PLL_C2S(VERSION_MINOR) "." \
        PLL_C2S(VERSION_PATCH)

/* MTOA maps matrix indexes to array index using row-major order */
/* N1=totNoInd*noLoci N2=totNoInd I=ind L=locus A=allele */

#define MTOA(I,L,A,N1,N2) ((N1)*(A))+((N2)*(L))+(I) 

extern char version[];  /* global variable: name of program using library -- for use with stderr */

struct indiv
{
  char indLabel[MAXLN];
  char popLabel[MAXLN];
  char locusLabel[MAXLN];
  char allele1[MAXALLNM];
  char allele2[MAXALLNM];
};

typedef struct data_params
{
  unsigned int noPops;
  unsigned int noLoci;
  int noAlleles[MAXLOCI][2];
  unsigned int noInd[MAXPOP];
  int totNoInd;
  char** popNames;
  char** locusNames;
  char** indNames[MAXPOP];
} datapar;

typedef struct data_hash
{
  GHashTable* popKeys;
  GHashTable* indKeys[MAXPOP];
  GHashTable* lociKeys;
  GHashTable* alleleKeys[MAXLOCI];
} dhash;

gboolean addKey(GHashTable* hash, char* mykey, int index);

int noKeys(GHashTable* hash);

int keyToIndex(GHashTable* hash, char* mykey);

void fillheader(const char version[]);

void show_header();

int cstring_cmp(const void *a, const void *b);

void prMemSz(unsigned int x);

int isMissing(char* x);

void get_line_checkdata(FILE* inputFile, struct indiv* ind, const regex_t regex);

void get_line(FILE* inputFile, struct indiv* ind);
  
void fillData(FILE* inputFile, int* dataArray, dhash* dh, datapar* dpar);

void readGData(FILE *inputFile, dhash* dh, datapar* dpar);
