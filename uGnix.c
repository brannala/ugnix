#ifndef STDLIBS
#define STDLIBS
#include <stdio.h>
#include <stdlib.h>
#include <gmodule.h> 
#include "uGnix.h"
#endif

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


void iterator(gpointer key, gpointer value, gpointer user_data)
{
  printf(user_data, key, GPOINTER_TO_INT(value));
}

gboolean addKey(GHashTable* hash, char* mykey, int index)
{
  gboolean y=FALSE;
  if(!g_hash_table_contains(hash, mykey))
    return g_hash_table_insert(hash, mykey, GINT_TO_POINTER(index));
  else
    return y;
}

int noKeys(GHashTable* hash)
{
  return g_hash_table_size(hash);
}

int keyToIndex(GHashTable* hash, char* mykey)
{
  return GPOINTER_TO_INT(g_hash_table_lookup(hash, mykey));
}

void printKeys(GHashTable* hash,char* phrase)
{
  g_hash_table_foreach(hash, (GHFunc)iterator, phrase);
}

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

void getIndNames( struct indiv* genoTypes, GHashTable* indKeys[], GHashTable* popKeys, gchar** popNames, int noPops, int* noInds )
{
  int totNoInds=0; // give each individual a unique integer index
  struct indiv* genos;
  for(int i=0; i<noPops; i++)
    {
      indKeys[i] = g_hash_table_new(g_str_hash, g_str_equal);
      genos = genoTypes;
      int currNoInds=0; // keep count of number of individuals in each population
      genos = genos->next;
      while(genos->next != NULL) 
	{
	  if(!strcmp(popNames[i],genos->popLabel))
	    {
	      if(addKey(indKeys[i],genos->indLabel,totNoInds))
		{ currNoInds++; totNoInds++; }
	    }
	  genos=genos->next;
	}
      noInds[i]=currNoInds;
    }
}

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

void getAlleleNames( struct indiv* genoTypes, GHashTable* alleleKeys[], gchar** locusNames, unsigned int noLoci, int** noAlleles )
{
  struct indiv* genos;
  for(int i=0; i<noLoci; i++)
    {
      alleleKeys[i] = g_hash_table_new(g_str_hash, g_str_equal);
      genos = genoTypes;
      int currNoAlleles=1; // keep count of number of alleles at each locus
      genos = genos->next;
      while(genos->next != NULL) 
	{
	  if(!strcmp(locusNames[i],genos->locusLabel))
	    {
	      if(addKey(alleleKeys[i],genos->allele1,isMissing(genos->allele1) ? 0 : currNoAlleles))
		if(!isMissing(genos->allele1)) currNoAlleles++; 
	      if(addKey(alleleKeys[i],genos->allele2,isMissing(genos->allele2) ? 0 : currNoAlleles))
		if(!isMissing(genos->allele2)) currNoAlleles++; 
	    }
	  genos=genos->next;
	}
      noAlleles[i][0]=currNoAlleles-1;
    }
}

int isMissing(char* x)
{
  if(strcmp("0",x)!=0 && strcmp("?",x)!=0)
    return 0;
  else
    return 1;
}
