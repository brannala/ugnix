
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
int opt_print_help = 0; /* print help information */

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


int main(int argc, char **argv)
{
  struct indiv* genoTypes;
  datapar dpar;
  dpar.noPops = 0;
  dpar.totNoInd = 0;
  dhash dh;
  dh.popKeys = g_hash_table_new(g_str_hash, g_str_equal);
  dh.lociKeys = g_hash_table_new(g_str_hash, g_str_equal);

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
  if(optind < argc)
    strcpy(fileName,argv[optind]);
  else
    { cmd_help(); return 1; }
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
	  if(opt_print_default)
	    {
	      printf("no_pop: %d\t tot_no_ind: %d\n\n",dpar.noPops,dpar.totNoInd);
	      for(int i = 0; i < dpar.noPops; i++)
		{
		  printf("PopID: %s\t",dpar.popNames[i]);
		  printf("no_ind: %d\t no_loci: %d\n",dpar.noInd[i],dpar.noLoci);
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
		printf("PopID: %s\n",dpar.popNames[i]);
		printKeys(dh.indKeys[i],"IndID: %s (%d)\n");
	      } 
	  if(opt_print_loci)
	    {
	      for(int i=0; i<dpar.noLoci; i++)
		printf("LocID: %s\t no_alleles: %d\t missing: %s\n",dpar.locusNames[i],
		       dpar.noAlleles[keyToIndex(dh.lociKeys,dpar.locusNames[i])][0]-1,
		       dpar.noAlleles[keyToIndex(dh.lociKeys,dpar.locusNames[i])][1] ? "Y" : "N");
	    } 
	}
      fclose(inputFile);
    }
}


