#ifndef STDLIBS
#define STDLIBS
#include "uGnix.h"
#endif

#define DEBUG 0
#define PROG_FUNCTION "heter"

char progheader[100];

void fillheader(const char version[])
{
  snprintf(progheader, 80,
           "%s %s %s",
           PROG_NAME, PROG_VERSION, version);
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/bpp\n");
  fprintf(stdout,"\n");
}

/* read genotype file in BA3 format */

struct indiv* readGFile(FILE* inputFile) 
{
  int nvars;
  struct indiv* head;
  struct indiv* current;
  struct indiv* new;
  head = NULL;
  char oneLine[MAXLN];
  head = malloc(sizeof(struct indiv));
  head->next = NULL;
  current = head;
  while((fgets(oneLine,sizeof(oneLine),inputFile)) != NULL)
    {
      new = malloc(sizeof(struct indiv));
      new->next = NULL;
      current->next = new;
      current = new;
      nvars = sscanf(oneLine,"%s %s %s %s %s",current->indLabel,current->popLabel,
		     current->locusLabel,current->allele1,current->allele2);
      if(nvars != 5) { printf("Error: incorrect number of entries at: %s",oneLine); return(NULL); }
    }
  return(head);
}

/* get index of array from matrix coordinates */

int matToArr(int ind, int locus, int allele, datapar dpar)
{
  return(dpar.totNoInd*dpar.noLoci*allele+dpar.totNoInd*locus+ind);
}

void fillData(struct indiv* genoTypes, int* dataArray, dhash* dh, datapar dpar)
{
  struct indiv* genos;
  genos = genoTypes;
  genos = genos->next;
  while(genos->next != NULL) 
    {
      dataArray[matToArr(keyToIndex(dh->indKeys[keyToIndex(dh->popKeys, genos->popLabel)],
				    genos->indLabel),keyToIndex(dh->lociKeys,genos->locusLabel),0,dpar)] =
	keyToIndex(dh->alleleKeys[keyToIndex(dh->lociKeys,genos->locusLabel)],genos->allele1);
      dataArray[matToArr(keyToIndex(dh->indKeys[keyToIndex(dh->popKeys, genos->popLabel)],
				    genos->indLabel),keyToIndex(dh->lociKeys,genos->locusLabel),1,dpar)] =
      keyToIndex(dh->alleleKeys[keyToIndex(dh->lociKeys,genos->locusLabel)],genos->allele2); 
      genos = genos->next;
    }
}


/* used by printKeys to iterate all keys and values in hash */

void iterator(gpointer key, gpointer value, gpointer user_data)
{
  printf(user_data, key, GPOINTER_TO_INT(value));
}

/* try to add key and value. return FALSE if key already exists */

gboolean addKey(GHashTable* hash, char* mykey, int index)
{
  gboolean y=FALSE;
  if(!g_hash_table_contains(hash, mykey))
    return g_hash_table_insert(hash, mykey, GINT_TO_POINTER(index));
  else
    return y;
}

/* get number of keys in hash */

int noKeys(GHashTable* hash)
{
  return g_hash_table_size(hash);
}

/* convert key to index value */

int keyToIndex(GHashTable* hash, char* mykey)
{
  return GPOINTER_TO_INT(g_hash_table_lookup(hash, mykey));
}

/* print all key value pairs in hash */

void printKeys(GHashTable* hash,char* phrase)
{
  g_hash_table_foreach(hash, (GHFunc)iterator, phrase);
}

/* returns 1 (true) if missing data exist and o (false) otherwise */

int isMissing(char* x)
{
  if(strcmp("0",x)!=0 && strcmp("?",x)!=0)
    return 0;
  else
    return 1;
}

/* create a hash of popKeys, set number of populations (noPops) and */
/* return array of population names/keys */

void getPopNames( struct indiv* genoTypes, dhash* dh, datapar* dpar)
{
  struct indiv* genos = genoTypes;
  int currNoPops=0;
  genos = genos->next;
  while(genos->next != NULL) 
    {
      if(addKey(dh->popKeys,genos->popLabel,currNoPops))
	currNoPops++;
      genos=genos->next;
    }
  dpar->popNames = (gchar **) g_hash_table_get_keys_as_array(dh->popKeys,&(dpar->noPops));
}

/* create a hash of names of all individuals from all populations assigning */
/* each a index in (0,noInds-1), set array of noInds per population */

void getIndNames( struct indiv* genoTypes, dhash* dh, datapar* dpar)
{
  dpar->totNoInd=0; // give each individual a unique integer index
  struct indiv* genos;
  for(int i=0; i<dpar->noPops; i++)
    {
      dh->indKeys[keyToIndex(dh->popKeys,dpar->popNames[i])] =
	g_hash_table_new(g_str_hash, g_str_equal);
      genos = genoTypes;
      int currNoInds=0; // keep count of number of individuals in each population
      genos = genos->next;
      while(genos->next != NULL) 
	{
	  if(!strcmp(dpar->popNames[i],genos->popLabel))
	    {
	      if(addKey(dh->indKeys[keyToIndex(dh->popKeys,dpar->popNames[i])],
			genos->indLabel,dpar->totNoInd))
		{ currNoInds++; dpar->totNoInd = dpar->totNoInd + 1; }
	    }
	  genos=genos->next;
	}
      dpar->noInd[keyToIndex(dh->popKeys,dpar->popNames[i])]=currNoInds;
    }
}

/* Create a hash of locus names assigning each an integer index in (0,noLoci-1). */
/* Set noLoci and return array of locus names/keys */

void getLociNames( struct indiv* genoTypes, dhash* dh, datapar* dpar)
{
  struct indiv* genos;
  genos = genoTypes;
  int currNoLoci=0;
  genos = genos->next;
  while(genos->next != NULL) 
    {
      if(addKey(dh->lociKeys,genos->locusLabel,currNoLoci))
	currNoLoci++;
      genos=genos->next;
    }
  dpar->locusNames = (gchar **) g_hash_table_get_keys_as_array(dh->lociKeys,&(dpar->noLoci));
}

/* Create a hash of allele names for each locus and indicate whether missing data exists. */
/* Return array of locus-specific allele counts (+1) and 0/1 for missing data */

void getAlleleNames( struct indiv* genoTypes, dhash* dh, datapar* dpar)
{
  struct indiv* genos;
  for(int i=0; i<MAXLOCI; i++)
    {
      dpar->noAlleles[i][0] = 1;
      dpar->noAlleles[i][1] = 0;
    }
  for(int i=0; i<dpar->noLoci; i++)
    dh->alleleKeys[i] = g_hash_table_new(g_str_hash, g_str_equal);
  genos = genoTypes;
  genos = genos->next;
  while(genos->next != NULL) 
    {
	  if(addKey(dh->alleleKeys[keyToIndex(dh->lociKeys, genos->locusLabel)],
		    genos->allele1, isMissing(genos->allele1) ? 0 :
		    dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0]))
	    {
	      if(!isMissing(genos->allele1))
		dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0] =
		  dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0]+1;
	      else
		dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][1]=1;
	    }
	  if(addKey(dh->alleleKeys[keyToIndex(dh->lociKeys, genos->locusLabel)],
		    genos->allele2,isMissing(genos->allele2) ? 0 :
		    dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0]))
	    {
	      if(!isMissing(genos->allele2))
		dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0] =
		  dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0]+1;
	      else
		dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][1]=1;
	    }
	  genos=genos->next;
    }
}

/* Get data parameters and create hashes */

void getDataParams(struct indiv* genoTypes, dhash* dh, datapar* dpar)
{
	  getPopNames(genoTypes,dh,dpar); 
	  getIndNames(genoTypes,dh,dpar); 
	  getLociNames(genoTypes,dh,dpar);
	  getAlleleNames(genoTypes,dh,dpar); 
}
