#include<uGnix.h>
#include "kinship_data.h"


/* options */
int opt_print_inbreeding = 0; /* print inbreeding coefficient */
int opt_print_kinship = 0; /* print kinship coefficient  */
int opt_print_relatedness = 0; /* print coefficient of relatedness */
int opt_print_help = 0; /* print help information */
int opt_verbose = 0; /* print diagnostic information in output file */

FILE *fp;
char filename[MAX_FILE_NAME];
char version[] = "kinship";

static void program_help()
{
  printf("Usage: kinship [OPTION]... FILE\n\n");
  printf("Calculate measures of relatedness from pedigrees: \n"
	 "-i inbreeding coefficients\n"
	 "-k pairwise kinship coefficients\n"
	 "-r pairwise relatedness coefficients\n\n");
}

static void print_msg()
{
  printf("Usage: kinship [OPTION]... FILE\n Try 'kinship -h' for more information.\n\n");
}


static void print_output(int **individual,int **father,int **mother,int n, char result )
{
  int z, i, j, M, F;
  float *phi;
  double *r;
  float *f;
  phi = malloc(n * n * sizeof(float)); /* allocating for kinship */
  r = malloc(n * n * sizeof(double)); /* allocating for relatedness */
  f = malloc(n * sizeof(float)); /* allocating for inbreeding */
  z = 0;

  /* kinship coefficient calculation*/
  for(j = 0; j < n; ++j)
    {
      if(*(*mother+j) == 0 && *(*father+j) == 0)
	{
	  z = z + 1;
	}
    }
  for(j = 1; j <= n ; ++j)
    {
      for(i = 1; i <= j; ++i){
	if(i==j)
	  {
	    if(i<=z)

	      {
		*(phi + (i-1)*n + (i-1)) = 0.5;
	      }
	    else
	      {
		*(phi + (i-1)*n + (i-1)) = 0.5 + 0.5 * *(phi + (*(*mother+(i-1))-1)*n + (*(*father+(i-1))-1)) ;
	      }
	  }
	else
	  {
	    if(i<=z && j<=z)
	      {
		*(phi + (i-1)*n + (j-1)) = 0 ;
		*(phi + (j-1)*n +(i-1)) = 0;
	      }
	    else
	      {
		*(phi + (i-1)*n + (j-1)) = 0.5* *(phi + (i-1)*n + (*(*mother+(j-1))-1)) +0.5* *(phi + (i-1)*n + (*(*father+(j-1))-1));
		*(phi + (j-1)*n +(i-1)) = *(phi + (i-1)*n +(j-1));
	      }
	  }
      }
    }

  /* inbreeding coefficient calculation */
  for(j = 0 ; j < n ; ++j)
    {
      if(*(*mother+j) == 0 && *(*father+j) == 0)
	{
	  *(f + j) = 0;
	}
      else
	{
	  M = *(*mother+j);
	  F = *(*father+j);
	  *(f + j) = *(phi + (M-1)*n + (F-1));
	}
    }

  /* coefficient of relatedness calculation */
  for(j = 0; j < n; ++j)
    {
      for(i = 0; i < n; ++i)
	{
	  *(r + i*n + j) = *(phi + i*n + j) * 2/sqrt((1 + *(f+i)) * (1 + *(f+j)));
	}
    }

  /* print kinship coefficient between a pair of individuals */
  if(result == 'k')
    {
      printf("Kinship Coefficients By Index\n");
      printf("----------------------------------->\n");
      printf("        ");
      for(j = 0; j < n; ++j)
	printf("%-8d ",j+1);
      printf("\n");
      for(j = 0; j < n; ++j)
	{
	  printf("%-8d",j+1);
	  for(i = 0; i < n; ++i)
	    {
	      printf("%-7f ", *(phi + i*n + j));
	    }
	  printf("\n");
	}
      printf("----------------------------------->\n");
      printf("\n");
    }
  
 /* print inbreeding coefficient for each individual */
  else if(result == 'i')
    {
      printf("Inbreeding Coefficients By Index\n");
      printf("------------------------------------>\n");
      for(j = 0; j < n ; ++j)
	{
	  printf("%-8d ",j+1);
	  printf("%-8f\n",*(f+j));
	}
      printf("----------------------------------->\n");
      printf("\n");
    }

  /* print coefficient of relatedness between a pair of individuals */
  else if(result == 'r')
    {
      printf("Relatedness Coefficients By Index\n");
      printf("----------------------------------->\n");
      printf("        ");
      for(j = 0; j < n; ++j)
	printf("%-8d ",j+1);
      printf("\n");
      for(j = 0; j < n; ++j)
	{
	  printf("%-8d",j+1);
	  for(i = 0; i < n; ++i)
	    {
	      printf("%-7f ", *(r + i*n + j));
	    }
	  printf("\n");
	}
      printf("----------------------------------->\n");
      printf("\n");
    }
  free(phi);
  free(r);
  free(f);
}



int main(int argc, char **argv)
{
  bool inputFromFile = false;
  int c;

  fillheader(version);
  show_header();

  while((c = getopt(argc, argv, "hikrv")) != -1)
    {
      switch(c)
	{
	case 'i':
	  opt_print_inbreeding = 1;
	  break;
	case 'k':
	  opt_print_kinship = 1;
	  break;
	case 'r':
	  opt_print_relatedness = 1;
	  break;
	case 'h':
	  opt_print_help = 1;
	  break;
	case 'v':
	  opt_verbose = 1;
	  break;
	case '?':
	  if (isprint (optopt))
	    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
	  else
	    fprintf(stderr,"Unknown option character `\\x%x'.\n",optopt);
	  return 1;
	}
    }


  if(optind < argc)
    {
	  strcpy(filename,argv[optind]);
	  fp = fopen(filename,"r+");
	  inputFromFile = true;
	  if (fp == NULL)
	    {
	      printf("%s: Could not open file %s\n",argv[0],filename);
	      return 1;
	    }
    }
  else
    {
      if(opt_print_help)
	program_help();
      else
	print_msg();
      return 1;
    }

  if(opt_print_help)
    {
      program_help();
      exit(0);
    }
  /* Indexing the individuals */
  int *E, *index;
  E = malloc(sizeof(int));
  index = malloc(sizeof(int));
  int n;
  int *individual, *father, *mother;
  family *head;

  /* All kinds of error check */
  format_error_check(fp,E); /* Checking the format of the data file*/
  printf("\n");
  rewind(fp);
  if(*E != 0) { fprintf(stderr,"The number of error(s):  %d\n",*E); exit(1);}
  head = Linked_List_Creation(fp); 
  logical_error_check(head); /* Checking all types of logical error in the data file*/

  /* Print the individuals with the correct index */
  if(*E == 0)
    {
      arrange(opt_verbose,head,&individual,&father,&mother,index); 
      n = (*index);
    }

  /* print inbreeding coefficient option */
  if(opt_print_inbreeding)
    {
      print_output(&individual,&father,&mother,n,'i');
    }

  /* print kinship coefficient option */
  if(opt_print_kinship)
    {
      print_output(&individual,&father,&mother,n,'k');
    }

  /* print relatedness coefficient option */  
  if(opt_print_relatedness)
    {
      print_output(&individual,&father,&mother,n,'r');
    }
  free(individual);
  free(father);
  free(mother);
  free(E);
  free(index);

  if(inputFromFile)
    {
      fclose(fp);
    }
  return 0;
}
