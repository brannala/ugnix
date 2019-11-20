#define MAXLN 100 
#define MAXIND 10000
#define MAXPOP 500
#define MAXLOCI 100000

struct indiv
{
  char indLabel[MAXLN];
  char popLabel[MAXLN];
  char locusLabel[MAXLN];
  char allele1[MAXLN];
  char allele2[MAXLN];
  struct indiv* next;
};
  
struct indiv* readGFile(FILE* inputFile);

int matToArr(int ind, int locus, int allele, int noInd, int noLoci);

void fillDataMatrix(struct indiv* genoTypes, int** dataMat, GHashTable* indKeys[], GHashTable* lociKeys, GHashTable* alleleKeys[], GHashTable* popKeys, int noLoci, int totNoInd);

void iterator(gpointer key, gpointer value, gpointer user_data);

gboolean addKey(GHashTable* hash, char* mykey, int index);

int noKeys(GHashTable* hash);

int keyToIndex(GHashTable* hash, char* mykey);

void printKeys(GHashTable* hash, char* phrase);

gchar** getPopNames( struct indiv* genoTypes, GHashTable* popKeys, unsigned int* noPops);

void getIndNames( struct indiv* genoTypes, GHashTable* indKeys[], GHashTable* popKeys, gchar** popNames, int noPops, int* noInd, int* totNoInds);

gchar** getLociNames( struct indiv* genoTypes, GHashTable* lociKeys, GHashTable* popKeys, gchar** popNames, int noPops, unsigned int* noLoci );

void getAlleleNames( struct indiv* genoTypes, GHashTable* alleleKeys[], GHashTable* locusKeys, unsigned int noLoci, int noAlleles[MAXLOCI][2] );

int isMissing(char* x);
