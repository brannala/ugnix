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
  fprintf(stdout, "https://github.com/ugnix\n");
  fprintf(stdout,"\n");
}

/* read one line of input file into struct indiv */
/* used for parsing by fillData() and readGdata */

void get_line(FILE* inputFile, struct indiv* ind)
{
  long int lineNo=1;
  int nvars;
  char oneLine[1000];
  if((fgets(oneLine,sizeof(oneLine),inputFile)) != NULL)
    {
      nvars = sscanf(oneLine,"%s %s %s %s %s",ind->indLabel,ind->popLabel,
		     ind->locusLabel,ind->allele1,ind->allele2);
      if(nvars != 5)
	{
	  fprintf(stderr,"\nError: At line %ld, %d != 5 is an incorrect number of entries at: %s\n",
		  lineNo,nvars,strtok(oneLine,"\n"));
	  exit(1);
	}
      ++lineNo;
    }
  else
    strcpy(ind->popLabel,"eof");
};

void fillData(FILE* inputFile, int* dataArray, dhash* dh, datapar* dpar)
{
  int n1 = dpar->totNoInd*dpar->noLoci;
  int n2 = dpar->totNoInd;
  struct indiv* genos;
  genos = malloc(sizeof(struct indiv));
  fprintf(stderr,"Generating data structures...\n\n");
  rewind(inputFile);
  get_line(inputFile,genos);
  while(strcmp(genos->popLabel,"eof"))
    {
      dataArray[MTOA(keyToIndex(dh->indKeys[keyToIndex(dh->popKeys, genos->popLabel)],
				genos->indLabel),keyToIndex(dh->lociKeys,genos->locusLabel),0,n1,n2)] =
	keyToIndex(dh->alleleKeys[keyToIndex(dh->lociKeys,genos->locusLabel)],genos->allele1);
      dataArray[MTOA(keyToIndex(dh->indKeys[keyToIndex(dh->popKeys, genos->popLabel)],
				genos->indLabel),keyToIndex(dh->lociKeys,genos->locusLabel),1,n1,n2)] =
      keyToIndex(dh->alleleKeys[keyToIndex(dh->lociKeys,genos->locusLabel)],genos->allele2); 
      get_line(inputFile,genos);
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
  if(strcmp("0",x)!=0 && strcmp("?",x)!=0 && strcmp(".",x)!=0)
    return 0;
  else
    return 1;
}

/* create a hash of popKeys, set number of populations (noPops) and */
/* return array of population names/keys */
/* create a hash of names of all individuals from all populations assigning */
/* each a index in (0,noInds-1), set array of noInds per population */
/* Create a hash of locus names assigning each an integer index in (0,noLoci-1). */
/* Set noLoci and return array of locus names/keys */
/* Create a hash of allele names for each locus and indicate whether missing data exists. */
/* Return array of locus-specific allele counts (+1) and 0/1 for missing data */
/* Get data parameters and create hashes */

void readGData(FILE* inputFile, dhash* dh, datapar* dpar)
{
  int currNoPops=0;
  int currNoLoci=0;
  struct indiv* genos;
  genos = malloc(sizeof(struct indiv));
  if(inputFile!=NULL)
    {
      fprintf(stderr,"Getting populations and locus labels...\n");
      get_line(inputFile,genos);      
      while(strcmp(genos->popLabel,"eof"))
	{
	  char* tmp_poplabel;
	  tmp_poplabel = strdup(genos->popLabel);
	  if(addKey(dh->popKeys,tmp_poplabel,currNoPops)) /* get population label */
	    currNoPops++;
	  else
	    free(tmp_poplabel);
	  char* tmp_locuslabel;
	  tmp_locuslabel = strdup(genos->locusLabel);
	  if(addKey(dh->lociKeys,tmp_locuslabel,currNoLoci)) /* get locus label */
	    currNoLoci++;
	  else
	    free(tmp_locuslabel);
	  get_line(inputFile,genos);
	}
      dpar->locusNames = (gchar **) g_hash_table_get_keys_as_array(dh->lociKeys,&(dpar->noLoci));
      dpar->popNames = (gchar **) g_hash_table_get_keys_as_array(dh->popKeys,&(dpar->noPops));
      fprintf(stderr,"Getting individuals and allele labels...\n\n");
      rewind(inputFile);
      for(int i=0; i<MAXLOCI; i++)
	{
	  dpar->noAlleles[i][0] = 1;
	  dpar->noAlleles[i][1] = 0;
	}
      for(int i=0; i<dpar->noLoci; i++)
	dh->alleleKeys[i] = g_hash_table_new(g_str_hash, g_str_equal);
      dpar->totNoInd=0; // give each individual a unique integer index
      for(int i=0; i<dpar->noPops; i++) /* get individual labels in each population */
	  dh->indKeys[i] =
	    g_hash_table_new(g_str_hash, g_str_equal);

      get_line(inputFile,genos);
      while(strcmp(genos->popLabel,"eof"))
	    {
	      char* tmp_indlabel;
	      tmp_indlabel = strdup(genos->indLabel); 
	      if(addKey(dh->indKeys[keyToIndex(dh->popKeys,genos->popLabel)],
			tmp_indlabel,dpar->totNoInd))
		{
		  dpar->noInd[keyToIndex(dh->popKeys,genos->popLabel)]=dpar->noInd[keyToIndex(dh->popKeys,genos->popLabel)]+1;
		  dpar->totNoInd = dpar->totNoInd + 1;
		}
	      else
		free(tmp_indlabel);
	      char* tmp_allele1;
	      tmp_allele1 = strdup(genos->allele1); 
	      if(addKey(dh->alleleKeys[keyToIndex(dh->lociKeys, genos->locusLabel)],
			tmp_allele1, isMissing(genos->allele1) ? 0 :
			dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0]))
		{
		  if(!isMissing(genos->allele1))
		    dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0] =
		      dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0]+1;
		  else
		    dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][1]=1;
		}
	      else
		free(tmp_allele1);
	      char* tmp_allele2;
	      tmp_allele2 = strdup(genos->allele2); 
	      if(addKey(dh->alleleKeys[keyToIndex(dh->lociKeys, genos->locusLabel)],
			tmp_allele2,isMissing(genos->allele2) ? 0 :
			dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0]))
		{
		  if(!isMissing(genos->allele2))
		    dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0] =
		      dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][0]+1;
		  else
		    dpar->noAlleles[keyToIndex(dh->lociKeys, genos->locusLabel)][1]=1;
		}
	      else
		free(tmp_allele2);
	      get_line(inputFile,genos);
	    }
    }
  free(genos);
}
