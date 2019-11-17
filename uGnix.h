#define MAXLN 100 
#define MAXIND 10000
#define MAXPOP 500
#define MAXLOCI 50000

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
void iterator(gpointer key, gpointer value, gpointer user_data);
gboolean addKey(GHashTable* hash, char* mykey, int index);
int noKeys(GHashTable* hash);
int keyToIndex(GHashTable* hash, char* mykey);
void printKeys(GHashTable* hash, char* phrase);
gchar** getPopNames( struct indiv* genoTypes, GHashTable* popKeys, unsigned int* noPops); /* creates hash to all unique population names and calculates noPops */
void getIndNames( struct indiv* genoTypes, GHashTable* indKeys[], GHashTable* popKeys, gchar** popNames, int noPops, int* noInd ); /* creates hash to all unique indIDs for each population */
void getLociNames( struct indiv* genoTypes, GHashTable* lociKeys[], GHashTable* popKeys, gchar** popNames, int noPops, long int* noLoci ); /* create hash to all unique loci names of each population */ 
