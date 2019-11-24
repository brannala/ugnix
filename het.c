
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
int opt_print_ind = 0; /* print labels of individuals */
int opt_print_no_pop = 0; /* print number of populations */
int opt_print_loci = 0; /* print details for each locus */
int opt_print_default = 0; /* print default summary (population names, noInd, noLoci) */

FILE* inputFile;
char fileName[100];
char version[] = "het";

static void heter(int* dataArray,dhash* dh,datapar dpar, char type)
{
  int n1 = dpar.totNoInd*dpar.noLoci;
  int n2 = dpar.totNoInd;
  if(type == 'i')
    {
      gchar** indList[MAXPOP];
      for(int i=0; i<dpar.noPops; i++)
	{
	  printf("Population: %s\n",dpar.popNames[i]);
	  unsigned int nInds;
	  indList[keyToIndex(dh->popKeys,dpar.popNames[i])] = (gchar **) g_hash_table_get_keys_as_array(dh->indKeys[keyToIndex(dh->popKeys,dpar.popNames[i])],&nInds);
	  for(int j=0; j<nInds; j++)
	    {
	      int indIndex = keyToIndex(dh->indKeys[keyToIndex(dh->popKeys,dpar.popNames[i])],indList[keyToIndex(dh->popKeys,dpar.popNames[i])][j]);
	      int total_genotypes = 0;
	      int total_hets = 0;
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
		}
  	      printf("Indiv: %s Heter: %f\n",indList[keyToIndex(dh->popKeys,dpar.popNames[i])][j],(total_hets+0.0)/total_genotypes);
	    }
	}
    }
}


int main(int argc, char **argv)
{
  struct indiv* genoTypes;
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
  while((c = getopt(argc, argv, "lpi")) != -1)
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
	  getDataParams(genoTypes,&dh,&dpar);
	  dataArray = malloc((dpar.noLoci*dpar.totNoInd*2+1) * sizeof(int));
	  if(opt_print_default)
	    {
	      

	      
	      int n1 = dpar.totNoInd*dpar.noLoci;
	      int n2 = dpar.totNoInd;
	      fillData(genoTypes,dataArray,&dh,dpar);
	      heter(dataArray,&dh,dpar,'i');
	      /*
	      printf("dataArray[0]: %d dataArray[dpar.noLoci*dpar.totNoInd*2]: %d",
		     dataArray[0],dataArray[dpar.noLoci*dpar.totNoInd*2]);
	      for(int i=0; i<dpar.totNoInd; i++)
		for(int j=0; j<dpar.noLoci; j++)
		  {
		    printf("Indiv %d, locus %d, allele1: %d, allele2: %d\n",
			   i,j,dataArray[MTOA(i,j,0,n1,n2)],
			   dataArray[MTOA(i,j,1,n1,n2)]);
		    
			   }*/
	      for(int i = 0; i < dpar.noPops; i++)
	      {

		// printf("PopID: %s\t",popNames[i]);
		// printf("PopIndex: %d\n",keyToIndex(popKeys, popNames[i]));
		// printf("NoInd: %d\t NoLoci: %d\n",noInd[i],noLoci);

	      }
	    }
	}
    }
}
