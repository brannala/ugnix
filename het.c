
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

#include "uGnix.h"

/* options */
int opt_heter_ind = 0; /* print average heterozygosity for each individual */
int opt_heter_loci = 0; /* print average heterozygosity for each locus */
int opt_help = 0; /* print help message */
int opt_default = 0; /* print help message */

FILE* inputFile;
char fileName[100];
char version[] = "het";

static void print_msg()
{
  printf("Usage: het [OPTION]... FILE\n Try 'het -h' for more information.\n");
}

static void print_help()
{
  printf("Heterozygosity calculations: \n"
	 "-i print average heterozygosity for each individual\n"
	 "-l print average heterozygosity for each locus\n"
	 "-a print heterozygosity for each individual at each locus\n");
}

static void heter(int* dataArray,dhash* dh,datapar dpar, char type)
{
  int n1 = dpar.totNoInd*dpar.noLoci;
  int n2 = dpar.totNoInd;
  if(type == 'i')  /* print avg heter across loci for each individual */
    {
      gchar** indList[MAXPOP];
      for(int i=0; i<dpar.noPops; i++)
	{
	  unsigned int nInds;
	  indList[keyToIndex(dh->popKeys,dpar.popNames[i])] =
	    (gchar **) g_hash_table_get_keys_as_array(dh->indKeys[keyToIndex(dh->popKeys,dpar.popNames[i])],&nInds);
	  for(int j=0; j<nInds; j++)
	    {
	      int indIndex = keyToIndex(dh->indKeys[keyToIndex(dh->popKeys,dpar.popNames[i])],
					indList[keyToIndex(dh->popKeys,dpar.popNames[i])][j]);
	      int total_genotypes = 0;
	      int total_hets = 0;
	      double h1 = 0;
	      double sd_h1 = 0;
	      for(int k=0; k<dpar.noLoci; k++)
		{
		  int a1 = dataArray[MTOA(indIndex,k,0,n1,n2)];
		  int a2 = dataArray[MTOA(indIndex,k,1,n1,n2)];
		  if((a1!=0)&&(a2!=0))
		    {
		      total_genotypes++;
		      if(a1 != a2)
			total_hets++;
		    }
		  h1 = (total_hets+0.0)/total_genotypes;
		  sd_h1 = sqrt(h1*(1-h1)/total_genotypes);
		}
	      if(total_genotypes>0)
		printf("PopID: %s\tIndID: %s\tH1: %.3f +/- %.3f\n",dpar.popNames[i],
		     indList[keyToIndex(dh->popKeys,dpar.popNames[i])][j],h1,(sd_h1*1.96));
	      else
		printf("PopID: %s\tIndID: %s\tH1: %s +/- %s",dpar.popNames[i],
		       indList[keyToIndex(dh->popKeys,dpar.popNames[i])][j],"M","M");
	    }
	}
    }
}


int main(int argc, char **argv)
{
  bool inputFromFile=false;
  datapar dpar;
  dpar.noPops = 0;
  dpar.totNoInd = 0;
  dhash dh;
  dh.popKeys = g_hash_table_new(g_str_hash, g_str_equal);
  dh.lociKeys = g_hash_table_new(g_str_hash, g_str_equal);
  int* dataArray;

  fillheader(version);
  show_header();
  
  opterr = 0;
  int c;
  while((c = getopt(argc, argv, "alih")) != -1)
    switch(c)
      {
      case 'i':
	opt_heter_ind = 1;
	break;
       case 'l':
	opt_heter_loci = 1;
	break;
      case 'h':
	opt_help = 1;
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
  
  if(optind == 1) opt_default=1;
  if(optind < argc) /* additional arguments present: filename? */
    {
    strcpy(fileName,argv[optind]);
    inputFile = fopen(fileName,"r");
    inputFromFile=true;
    if( inputFile == NULL )
      {
	printf("%s: stat of %s failed: no such file\n",argv[0],fileName);
	return 1;
      }
    }
  else
    {
      if(isatty(STDIN_FILENO))  /* is input from a tty rather than a redirect or pipe? */
	{                       /* if from tty exit with help msg */
	  if(opt_help)
	    print_help();
	  else
	    print_msg();
	  return 1;
	}
    }
  if(inputFromFile)
    readGData(inputFile,&dh,&dpar);
  dataArray = malloc((dpar.noLoci*dpar.totNoInd*2+1) * sizeof(int));
  if(inputFromFile)
    fillData(inputFile,dataArray,&dh,&dpar);
  if(opt_heter_ind)
    {
      heter(dataArray,&dh,dpar,'i');
    }
  if(opt_heter_loci)
    {
      heter(dataArray,&dh,dpar,'l');
    }
  if(opt_help)
    {
      print_help();
    }
}

