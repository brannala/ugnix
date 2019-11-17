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

  struct indiv* genos;
  for(int i=0; i<noPops; i++)
    {
      indKeys[i] = g_hash_table_new(g_str_hash, g_str_equal);
      genos = genoTypes;
      int currNoInds=0;
      genos = genos->next;
      while(genos->next != NULL) 
	{
	  if(!strcmp(popNames[i],genos->popLabel))
	    {
	      if(addKey(indKeys[i],genos->indLabel,currNoInds))
		currNoInds++;
	    }
	  genos=genos->next;
	}
      noInds[i]=currNoInds;
    }
}

void getLociNames( struct indiv* genoTypes, GHashTable* lociKeys[], GHashTable* popKeys, gchar** popNames, int noPops, long int* noLoci )
{
  struct indiv* genos;
  for(int i=0; i<noPops; i++)
    {
      lociKeys[i] = g_hash_table_new(g_str_hash, g_str_equal);
      genos = genoTypes;
      int currNoLoci=0;
      genos = genos->next;
      while(genos->next != NULL) 
	{
	  if(!strcmp(popNames[i],genos->popLabel))
	    {
	      if(addKey(lociKeys[i],genos->locusLabel,currNoLoci))
		currNoLoci++;
	    }
	  genos=genos->next;
	}
      noLoci[i]=currNoLoci;
    }
}
