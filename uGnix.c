#ifndef STDLIBS
#define STDLIBS
#include <stdio.h>
#include <stdlib.h>
#include <gmodule.h> 
#include "uGnix.h"
#endif

/* read genotype file in BA3 format */
struct indiv* readGFile(FILE* inputFile) 
{
  int nvars;
  struct indiv* head;
  struct indiv* current;
  struct indiv* new;
  head = NULL;
  char oneLine[1000];
  head = malloc(sizeof(struct indiv));
  head->next = NULL;
  current = head;
  while((fgets(oneLine,sizeof(oneLine),inputFile)) != NULL)
    {
      new = malloc(sizeof(struct indiv));
      new->next = NULL;
      current->next = new;
      current = new;
      nvars = sscanf(oneLine,"%s %s %s %s %s",current->indLabel,current->popLabel,current->locusLabel,current->allele1,current->allele2);
      if(nvars != 5) { printf("Error: incorrect number of entries at: %s",oneLine); return(NULL); }
    }
  return(head);
}

/* get index of array from matrix coordinates */
int matToArr(int ind, int locus, int allele, int noInd, int noLoci)
{
  return(noInd*noLoci*allele+noInd*locus+ind);
}

void fillDataMatrix(struct indiv* genoTypes, int** dataMat, GHashTable* indKeys[], GHashTable* lociKeys, GHashTable* alleleKeys[], GHashTable* popKeys, int noLoci, int totNoInd)
{
  struct indiv* genos;
  genos = genoTypes;
  genos = genos->next;
  while(genos->next != NULL) 
    {
      printf("indiv: %s, locus: %s, allele1: %s, allele2: %s, arrInd: %d\n",genos->indLabel,genos->locusLabel,genos->allele1,genos->allele2,matToArr(keyToIndex(indKeys[keyToIndex(popKeys, genos->popLabel)], genos->indLabel),keyToIndex(lociKeys,genos->locusLabel),0,totNoInd,noLoci));
            printf("indiv: %s, locus: %s, allele1: %s, allele2: %s, arrInd: %d\n",genos->indLabel,genos->locusLabel,genos->allele1,genos->allele2,matToArr(keyToIndex(indKeys[keyToIndex(popKeys, genos->popLabel)], genos->indLabel),keyToIndex(lociKeys,genos->locusLabel),1,totNoInd,noLoci));
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

/* create a hash of popKeys, set number of populations (noPops) and return array of population names/keys */
gchar** getPopNames( struct indiv* genoTypes, GHashTable* popKeys, unsigned int* noPops)
{
  struct indiv* genos = genoTypes;
  int currNoPops=0;
  genos = genos->next;
  while(genos->next != NULL) 
    {
      if(addKey(popKeys,genos->popLabel,currNoPops))
	currNoPops++;
      genos=genos->next;
    }
  return (gchar **) g_hash_table_get_keys_as_array(popKeys,noPops);
}

/* create a hash of names of all individuals from all populations assigning each a index in (0,noInds-1), set array of noInds per population */
void getIndNames( struct indiv* genoTypes, GHashTable* indKeys[], GHashTable* popKeys, gchar** popNames, int noPops, int* noInds, int* totNoInds)
{
  *totNoInds=0; // give each individual a unique integer index
  struct indiv* genos;
  for(int i=0; i<noPops; i++)
    {
      indKeys[keyToIndex(popKeys,popNames[i])] = g_hash_table_new(g_str_hash, g_str_equal);
      genos = genoTypes;
      int currNoInds=0; // keep count of number of individuals in each population
      genos = genos->next;
      while(genos->next != NULL) 
	{
	  if(!strcmp(popNames[i],genos->popLabel))
	    {
	      if(addKey(indKeys[keyToIndex(popKeys,popNames[i])],genos->indLabel,*totNoInds))
		{ currNoInds++; *totNoInds = *totNoInds + 1; }
	    }
	  genos=genos->next;
	}
      noInds[keyToIndex(popKeys,popNames[i])]=currNoInds;
    }
}

/* create a hash of locus names assigning each an integer index in (0,noLoci-1). Set noLoci and return array of locus names/keys */
gchar** getLociNames( struct indiv* genoTypes, GHashTable* lociKeys, GHashTable* popKeys, gchar** popNames, int noPops, unsigned int* noLoci )
{
  struct indiv* genos;
  genos = genoTypes;
  int currNoLoci=0;
  genos = genos->next;
  while(genos->next != NULL) 
    {
      if(addKey(lociKeys,genos->locusLabel,currNoLoci))
	currNoLoci++;
      genos=genos->next;
    }
  return (gchar **) g_hash_table_get_keys_as_array(lociKeys,noLoci);
}

/* create a hash of allele names for each locus and indicate whether missing data exists. Return array of locus-specific allele counts (+1) and 0/1 for missing data */
void getAlleleNames( struct indiv* genoTypes, GHashTable* alleleKeys[], GHashTable* locusKeys, unsigned int noLoci, int noAlleles[MAXLOCI][2] )
{
  struct indiv* genos;
  for(int i=0; i<MAXLOCI; i++)
    {
      noAlleles[i][0] = 1;
      noAlleles[i][1] = 0;
    }
  for(int i=0; i<noLoci; i++)
    alleleKeys[i] = g_hash_table_new(g_str_hash, g_str_equal);
  genos = genoTypes;
  genos = genos->next;
  while(genos->next != NULL) 
    {
	  if(addKey(alleleKeys[keyToIndex(locusKeys, genos->locusLabel)],genos->allele1,isMissing(genos->allele1) ? 0 : noAlleles[keyToIndex(locusKeys, genos->locusLabel)][0]))
	    {
	      if(!isMissing(genos->allele1))
		noAlleles[keyToIndex(locusKeys, genos->locusLabel)][0]=noAlleles[keyToIndex(locusKeys, genos->locusLabel)][0]+1;
	      else
		noAlleles[keyToIndex(locusKeys, genos->locusLabel)][1]=1;
	    }
	  if(addKey(alleleKeys[keyToIndex(locusKeys, genos->locusLabel)],genos->allele2,isMissing(genos->allele2) ? 0 : noAlleles[keyToIndex(locusKeys, genos->locusLabel)][0]))
	    {
	      if(!isMissing(genos->allele2))
		noAlleles[keyToIndex(locusKeys, genos->locusLabel)][0]=noAlleles[keyToIndex(locusKeys, genos->locusLabel)][0]+1;
	      else
		noAlleles[keyToIndex(locusKeys, genos->locusLabel)][1]=1;
	    }
	  genos=genos->next;
    }
}


