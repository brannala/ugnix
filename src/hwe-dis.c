
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
#include "uGnix.h"
#endif

/* options */
int opt_heter_ind = 0; /* print average heterozygosity for each individual */
int opt_heter_loci = 0; /* print average heterozygosity for each locus */
int opt_help = 0; /* print help message */
int opt_default = 0; /* print help message */

FILE* inputFile;
char fileName[100];
char version[] = "hwe-dis";

static void print_msg()
{
  printf("Usage: hwe-dis [OPTION]... FILE\n Try 'hwe-dis -h' for more information.\n");
}

static void print_help()
{
  printf("Additive Disequilibrium Test for Departure from Hardy-Weinberg Equilibrium: \n"
	 "-h print help\n");
}

static void add_diseq(int* dataArray, dhash* dh, datapar dpar)
{
  int n1 = dpar.totNoInd*dpar.noLoci;
  int n2 = dpar.totNoInd;
      for(int i=0; i<dpar.noPops; i++)
	{
	  int popIndx = keyToIndex(dh->popKeys,dpar.popNames[i]);
	  for(int j=0; j<dpar.noLoci; j++)
	    {
	      int locIndx = keyToIndex(dh->lociKeys,dpar.locusNames[j]);
	      int total_genotypes = 0;
	      int total_hets = 0;
	      double h1 = 0;
	      double sd_h1 = 0;
	      int targetAllele=-100;
	      int n11=0;
	      int n12=0;
	      double diseq=0.0;
	      double alleleFreq=0.0;
	      double sdev=0.0;
	      int significant=0;
	      int too_few_counts=0;
	      for(int k=0; k<dpar.noInd[popIndx]; k++)
		{
		  int indIndx = keyToIndex(dh->indKeys[popIndx],dpar.indNames[popIndx][k]); 
		  int a1 = dataArray[MTOA(indIndx,locIndx,0,n1,n2)];
		  int a2 = dataArray[MTOA(indIndx,locIndx,1,n1,n2)];
		  if(targetAllele==-100)
		    targetAllele=a1;
		  if((a1!=0)&&(a2!=0))
		    {
		      total_genotypes++;
		      if((targetAllele==a1)&&(targetAllele==a2))
			n11++;
		      else
			if((targetAllele==a1)||(targetAllele==a2))
			  n12++;
		    }
		}
	      if((n11 < 5)||(total_genotypes-n11 < 5))
		too_few_counts=1;
	      diseq = n11/(total_genotypes+0.0) - ((2.0*n11+n12)/(2.0*total_genotypes))*((2.0*n11+n12)/(2.0*total_genotypes));
	      alleleFreq = (2.0*n11+n12)/(2.0*total_genotypes);
	      sdev = sqrt(alleleFreq*alleleFreq*(1-alleleFreq)*(1-alleleFreq)/total_genotypes);
	      if(diseq >= 0)
		{
		  if((diseq - 1.96*sdev) < 0.0)
		    significant = 0;
		  else
		    significant = 1;
		}
	      else
		{
		  if((diseq + 1.96*sdev) > 0.0)
		    significant = 0;
		  else
		    significant = 1;
		}
	      if((total_genotypes>0)&&(!too_few_counts))
		{
		  printf("PopID: %s\tlocID: %s\tD: %.5f +/- %.5f ",dpar.popNames[i],
		       dpar.locusNames[j],diseq,1.96*sdev);
		  if(significant)
		    printf("*\n");
		  else
		    printf("\n");
		}
	      else
		printf("PopID: %s\tlocID: %s\tD: %s\n",dpar.popNames[i],
		       dpar.locusNames[j],"M");
	    }
	}
}

int main(int argc, char **argv)
{
  bool inputFromFile=false;
  datapar dpar = {.noPops=0, .totNoInd=0};
  dhash dh = {.popKeys = g_hash_table_new(g_str_hash, g_str_equal),
	      .lociKeys = g_hash_table_new(g_str_hash, g_str_equal)}; 
  int* dataArray;
  fillheader(version);
  show_header();
  opterr = 0;
  int c;
  while((c = getopt(argc, argv, "h")) != -1)
    switch(c)
      {
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
  unsigned int datasize = (dpar.noLoci*dpar.totNoInd*2+1); 
  if((dataArray = calloc(datasize,sizeof(int)))==NULL)
    {
      fprintf(stderr,"%s: out of memory! exiting gracefully...\n",version);
      exit(1);
    }
  prMemSz(datasize*sizeof(int));
  fillData(inputFile,dataArray,&dh,&dpar);

  if(opt_help)
    {
      print_help();
    }
  add_diseq(dataArray,&dh,dpar);

  /* cleanup memory */
  free(dataArray);
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

