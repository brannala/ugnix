#ifndef STDLIBS
#define STDLIBS
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <gmodule.h>
#include "uGnix.h"
#include <unistd.h>
#endif

int main(int argc, char **argv)
{
  int printInd = 0;
  int printNoPop = 0;
  int printLoci = 0;
  int printDefault = 0;
  int c;
  unsigned int noPops = 0;
  int noInd[MAXPOP];
  int totNoInd;
  unsigned int noLoci;
  char** popNames;
  char** locusNames;
  struct indiv* genoTypes;
  char fileName[100];
  int** dataMat;
  GHashTable* popKeys;
  popKeys = g_hash_table_new(g_str_hash, g_str_equal);
  GHashTable* indKeys[MAXPOP];
  GHashTable* lociKeys;
  lociKeys = g_hash_table_new(g_str_hash, g_str_equal);
  GHashTable* alleleKeys[10000];
  int noAlleles[MAXLOCI][2];
  FILE* inputFile;
  
  opterr = 0;

  while((c = getopt(argc, argv, "lpi")) != -1)
    switch(c)
      {
      case 'i':
	printInd = 1;
	break;
      case 'p':
	printNoPop = 1;
	break;
      case 'l':
	printLoci = 1;
	break;
      case '?':
        if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr,"Unknown option character `\\x%x'.\n",optopt);
        return 1;
      default:
	abort();
      }
  
  if(optind == 1) printDefault=1;
  if(optind < argc)
    strcpy(fileName,argv[optind]);
  else
    { printf("Missing filename argument!\n"); return 1; }
  inputFile = fopen(fileName,"r");
  if( inputFile == NULL )
    printf("%s: stat of %s failed: no such file\n",argv[0],fileName);
  else
    {
      genoTypes = readGFile(inputFile);
      if(genoTypes == NULL)
	  return 1;
      else
	{
	  popNames = getPopNames(genoTypes,popKeys,&noPops); 
	  getIndNames(genoTypes,indKeys,popKeys,popNames,noPops,noInd,&totNoInd); 
	  locusNames = getLociNames(genoTypes,lociKeys,popKeys,popNames,noPops,&noLoci);
	  getAlleleNames(genoTypes,alleleKeys,lociKeys,noLoci,noAlleles);
	  if(printDefault)
	    for(int i = 0; i < noPops; i++)
	      {
		fillDataMatrix(genoTypes,dataMat,indKeys,lociKeys,alleleKeys,popKeys,noLoci,totNoInd);
		// printf("PopID: %s\t",popNames[i]);
		// printf("PopIndex: %d\n",keyToIndex(popKeys, popNames[i]));
		// printf("NoInd: %d\t NoLoci: %d\n",noInd[i],noLoci);

	      }
	}
    }
}
