
/*
    Copyright (C) 2019 Bruce Rannala

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#endif

/* options */
int opt_print_ind = 0;     /* print labels of individuals */
int opt_print_no_pop = 0;  /* print number of populations */
int opt_print_loci = 0;    /* print details for each locus */
int opt_print_default = 0; /* print default summary (population names, noInd, noLoci) */
int opt_print_help = 0;    /* print help information */

FILE* inputFile;
char fileName[100];
char version[] = "gsum";

static void cmd_help()
{
  /*       0         1         2         3         4         5         6         7          */
  /*       01234567890123456789012345678901234567890123456789012345678901234567890123456789 */

  fprintf(stderr,
          "Usage: %s [OPTIONS]... FILE \n"
	  "List information about the FILE (basic summary of data by default).", version);
  fprintf(stderr,
          "\n"
          "General options:\n"
          "  -h                 display help information\n"
          "  -v                 display version information\n"
          "  -i                 summarize individuals\n"
          "  -p                 summarize populations\n"
          "  -l                 summarize loci\n"
	  "Notice: FILE must be in BA3/Immanc format.\n"
          "\n"
         );

  /*       0         1         2         3         4         5         6         7          */
  /*       01234567890123456789012345678901234567890123456789012345678901234567890123456789 */
}

/* prints formatted alleleIDs in alphabetical order */

static void printSortedAlleles(GHashTable* hash,char* phrase)
{
  unsigned int len;
  char** keyArray = (gchar **) g_hash_table_get_keys_as_array(hash,&len);
  qsort(keyArray,len,sizeof(char *),cstring_cmp);
  printf("%s",phrase);
  for(unsigned int i=0; i<len; i++)
    {
      printf(" %s",keyArray[i]);
    }
  g_free(keyArray);
}

int main(int argc, char **argv)
{
  bool inputFromFile=false;
  datapar dpar = {.noPops=0, .totNoInd=0};
  dhash dh = {.popKeys = g_hash_table_new(g_str_hash, g_str_equal),
	      .lociKeys = g_hash_table_new(g_str_hash, g_str_equal)}; 
  fillheader(version);
  show_header();
  opterr = 0;
  int c;
  while((c = getopt(argc, argv, "lpih")) != -1)
    switch(c)
      {
      case 'i':
	opt_print_ind = 1;
	break;
      case 'p':
	opt_print_no_pop = 1;
	break;
      case 'l':
	opt_print_loci = 1;
	break;
      case 'h':
	opt_print_help = 1;
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
  if(optind == 1) opt_print_default=1;
  if(optind < argc)  /* additional arguments present: filename? */
    {
    strcpy(fileName,argv[optind]);
    inputFile = fopen(fileName,"r");
    inputFromFile=true;
    if( inputFile == NULL )
      {
	printf("%s: stat of %s failed: no such file\n",version,fileName);
	exit(1);
      }
    }
  else
    {
      if(isatty(STDIN_FILENO))  /* is input from a tty rather than a redirect or pipe? */
	{                       /* if from tty exit with help msg */
	  cmd_help();
	  return 1;
	}
    }
  if(inputFromFile)
    readGData(inputFile,&dh,&dpar);
  else
    ; /* readGData(stdin,&dh,&dpar); broken. need to buffer to file because we can't rewind stdin. */
  if(opt_print_default)
    {
      printf("no_pop: %d\t tot_no_ind: %d\n\n",dpar.noPops,dpar.totNoInd);
      for(int i = 0; i < dpar.noPops; i++)
	{
	  printf("PopID: %s\t",dpar.popNames[i]);
	  printf("no_ind: %d\t no_loci: %d\n",dpar.noInd[keyToIndex(dh.popKeys,dpar.popNames[i])],dpar.noLoci);
	}
      printf("\n");
    }
  if(opt_print_help)
    {
      cmd_help();
      return 1;
    }
  if(opt_print_no_pop)
    for(int i = 0; i < dpar.noPops; i++)
      {
	printf("PopID: %s\n",dpar.popNames[i]);
      }
  if(opt_print_ind)
    for(int i = 0; i < dpar.noPops; i++)
      {
	int popIndx = keyToIndex(dh.popKeys,dpar.popNames[i]);
	char formatStr[50];
	strcpy(formatStr,"PopID: ");
	strcat(formatStr,dpar.popNames[i]);
	strcat(formatStr,"\tIndID: ");
	for(int i=0; i<dpar.noInd[popIndx]; i++)
	  printf("%s %s\n",formatStr,dpar.indNames[popIndx][i]);
      } 
  if(opt_print_loci)
    {
      for(int i=0; i<dpar.noLoci; i++)
	{
	  printf("LocID: %s\t no_alleles: %d\t missing: %s",dpar.locusNames[i],
		 dpar.noAlleles[keyToIndex(dh.lociKeys,dpar.locusNames[i])][0]-1,
		 dpar.noAlleles[keyToIndex(dh.lociKeys,dpar.locusNames[i])][1] ? "Y" : "N");
	  int currNoAlleles = dpar.noAlleles[keyToIndex(dh.lociKeys,dpar.locusNames[i])][0]-1;
	  if( currNoAlleles > 0)
	    {
	      printSortedAlleles(dh.alleleKeys[keyToIndex(dh.lociKeys,dpar.locusNames[i])],"\t alleles: ");
	      printf("\n");
	    }
	  else
	    printf("\n");
	}
    } 
  /* cleanup memory */
  g_hash_table_destroy(dh.popKeys);
  g_hash_table_destroy(dh.lociKeys);
  for(int i=0; i<dpar.noPops; i++)
    g_hash_table_destroy(dh.indKeys[i]);
  for(int i=0; i<dpar.noLoci; i++)
    g_hash_table_destroy(dh.alleleKeys[i]);
  if(inputFromFile)
    fclose(inputFile);

  return 0;
}



