#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <gmodule.h>
#include <unistd.h>

#define MAXLN 1000 
#define MAXIND 10000
#define MAXPOP 500
#define MAXLOCI 100000

/* constants */

#define PROG_NAME "uGnix"

#define PLL_STRING(x) #x
#define PLL_C2S(x) PLL_STRING(x)

#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_PATCH 1

#define PROG_VERSION "v" PLL_C2S(VERSION_MAJOR) "." PLL_C2S(VERSION_MINOR) "." \
        PLL_C2S(VERSION_PATCH)




struct indiv
{
  char indLabel[MAXLN];
  char popLabel[MAXLN];
  char locusLabel[MAXLN];
  char allele1[MAXLN];
  char allele2[MAXLN];
  struct indiv* next;
};

typedef struct data_params
{
  unsigned int noPops;
  unsigned int noLoci;
  int noAlleles[MAXLOCI][2];
  int noInd[MAXPOP];
  int totNoInd;
  char** popNames;
  char** locusNames;
} datapar;

typedef struct data_hash
{
  GHashTable* popKeys;
  GHashTable* indKeys[MAXPOP];
  GHashTable* lociKeys;
  GHashTable* alleleKeys[MAXLOCI];
} dhash;

struct indiv* readGFile(FILE* inputFile);

void fillheader(const char version[]);

void show_header();

int matToArr(int ind, int locus, int allele, datapar dpar);

void fillData(struct indiv* genoTypes, int* dataArray, dhash* dh, datapar dpar);

void iterator(gpointer key, gpointer value, gpointer user_data);

gboolean addKey(GHashTable* hash, char* mykey, int index);

int noKeys(GHashTable* hash);

int keyToIndex(GHashTable* hash, char* mykey);

void printKeys(GHashTable* hash, char* phrase);

void getPopNames( struct indiv* genoTypes, dhash* dh, datapar* dpar);

void getIndNames( struct indiv* genoTypes, dhash* dh, datapar* dpar);

void getLociNames( struct indiv* genoTypes, dhash* dh, datapar* dpar);

void getAlleleNames( struct indiv* genoTypes, dhash* dh, datapar* dpar);

void getDataParams(struct indiv* genoTypes, dhash* dh, datapar* dpar);

int isMissing(char* x);
