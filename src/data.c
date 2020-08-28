#include "kinship_data.h"
#include<assert.h>
#define BUFFER_SIZE   1000

int never_visited = 1;

family* makeFamily(char* indv, char* dad, char* mom)
{
  family* newFamily;
  if((newFamily = malloc(sizeof(family)))==NULL)
    { fprintf(stderr,"Oops, out of memory!"); exit(1); }

  strcpy(newFamily->indv,indv);
  strcpy(newFamily->dad,dad);
  strcpy(newFamily->mom,mom);
  newFamily->next = NULL;
  return(newFamily);
}

void iterator(gpointer key, gpointer value, gpointer user_data) {
  printf(user_data,*(gint*)value, key);
}

void format_error_check( FILE* f, int *e )
  {
    char *pattern;
    pattern = "^[[:blank:]]*([[:alnum:]_]{,20}[[:blank:]]+){2}([[:alnum:]_]{,20})([[:blank:]]*\n*)$|(^[[:blank:]]*\n$)";
    int t;
    regex_t re;
    char    buffer[BUFFER_SIZE];
    char row[BUFFER_SIZE];
    const char s[] = " \t";
    char* tok;
    int lineNo = 0;
    int error_count = 0;

    if ((t=regcomp( &re, pattern, REG_EXTENDED)) != 0) {
      regerror(t, &re, buffer, sizeof(buffer));
      fprintf(stderr,"%s (%s)\n",buffer,pattern);
      return;
    }
    while( fgets( buffer, BUFFER_SIZE, f ) != NULL ) {
      lineNo++;
      if( regexec( &re, buffer, 0, NULL, 0 ) != 0 ) {
	strcpy(row,buffer);
	int ncol = 0;
	tok = strtok(row,s);
	while (tok != 0)
	  {
	    ncol++;
	    tok = strtok(0,s);
	  }
	fprintf(stderr,"ERROR at Line %d: %s \t Incorrect number of entries: %d observed, but 3 columns expected\n",lineNo,buffer,ncol);
	error_count++;
      }
    }
    *e = error_count;
    regfree( &re );
  }


family* Linked_List_Creation(FILE *filename)
{
  char *pattern;
  pattern = "^[[:blank:]]*\n$";
  int t;
  regex_t re;
  char row[1000];
  family* head;
  family* current;
  char* i1;
  char* d1;
  char* m1;

  i1 = malloc(sizeof(char)*MAX_NAME);
  d1 = malloc(sizeof(char)*MAX_NAME);
  m1 = malloc(sizeof(char)*MAX_NAME);
  assert((t = regcomp(&re, pattern, REG_EXTENDED)) == 0);
  fgets(row, sizeof(row), filename);
  while(regexec(&re, row, 0, NULL, 0) == 0)
    {
      fgets(row, sizeof(row), filename);
    }
  sscanf(row,"%s %s %s", i1, d1, m1);
  head = makeFamily(i1,d1,m1);

  
  int first = 1;     
  int endoffile = 0;
  while ((fgets(row, sizeof(row), filename) != NULL) && !endoffile)
    {
      while(regexec(&re, row, 0, NULL, 0) == 0)
	{
	  if(fgets(row, sizeof(row), filename) == NULL)
	    {
	      endoffile = 1;
	      break;
	    }
	}
      if(!endoffile)
	{
	  sscanf(row,"%s %s %s", i1, d1, m1);
	  family* nextFamily = makeFamily(i1,d1,m1);

	  if (first == 1)
	    {
	      current = nextFamily;
	      head->next = current;
	      first = 0;
	    }
	  else
	    {
	      current->next = nextFamily;
	      current = current->next;
	    }
	}		
    }
  free(i1);
  free(d1);
  free(m1);
  regfree(&re);
  return head;	
}

void Hash_Function_creation(family *head, GHashTable* hash, GHashTable* hash_founder, int *n_founder, int *n_family_index, int node_count )
{
  family *current;
  current = head;
  int count;
  gint *Index;
  int current_Indv_Index = 0;
  while(current != NULL) /*********  FIRST PASS THROUGH LINKED LIST **********/
     {
       /*************** Indexing the founders first *********************/
       if(strcmp(current->mom,"0") == 0  && strcmp(current->dad,"0" ) == 0)
	{
	  current_Indv_Index++;
	  count = current_Indv_Index;
	  Index = g_new(gint,1);
	  *Index = count;
	  g_hash_table_insert(hash, current->indv, Index);
	  g_hash_table_insert(hash_founder, current->indv, Index);
	}
       current = current->next;
       
     }
  *n_founder = current_Indv_Index;
  current = head;
  while(current != NULL) /************************  SECOND PASS THROUGH LINKED LIST *********************/
    {
      /************* Indexing individuals both of whose parents are founders *****************/
      if(g_hash_table_contains(hash,current->indv) == FALSE && g_hash_table_contains(hash_founder,current->mom) == TRUE &&
	 g_hash_table_contains(hash_founder,current->dad) == TRUE )
	{
	  current_Indv_Index++;
	  count = current_Indv_Index;
	  Index = g_new(gint,1);
	  *Index = count;
	  g_hash_table_insert(hash, current->indv, Index);
	}
      current = current->next;
    }
  
  current = head;
  while(current != NULL) /************************  THIRD PASS THROUGH LINKED LIST *********************/
    {
      /* Indexing individuals when one of the parent is founder and the other nonfounder parent has been indexed already */
      if((g_hash_table_contains(hash,current->indv) == FALSE && g_hash_table_contains(hash_founder,current->mom) == TRUE &&
	  g_hash_table_contains(hash,current->dad) == TRUE) ||
	 (g_hash_table_contains(hash,current->indv) == FALSE && g_hash_table_contains(hash,current->mom) == TRUE &&
	  g_hash_table_contains(hash_founder,current->dad) == TRUE))
	{
	  current_Indv_Index++;
	  count = current_Indv_Index;
	  Index = g_new(gint,1);
	  *Index = count;
	  g_hash_table_insert(hash, current->indv, Index);
	}
      current = current->next;
    }
  
  while(current_Indv_Index < node_count)
    {
      current = head;
      while(current != NULL) /************************  FOURTH (and Last) PASS THROUGH LINKED LIST ******************** */
	{
	  /* When both parents have been indexed */
	  if(g_hash_table_contains(hash,current->indv) == FALSE && g_hash_table_contains(hash,current->mom) == TRUE &&
	     g_hash_table_contains(hash,current->dad) == TRUE)
	    {
	      current_Indv_Index++;
	      count = current_Indv_Index;
	      Index = g_new(gint,1);
	      *Index = count;
	      g_hash_table_insert(hash, current->indv, Index);
	    }
	  current = current->next;
	}      
    }
  *n_family_index = current_Indv_Index;
}

int node_count(family *head)
{
  int count = 0;
  while(head != NULL)
    {
      count++;
      head = head->next;
    }
  return(count);
}

void duplicate_row_check(family *head, GHashTable *hash)
{
  family *current;
  int count = 1;
  gint *Index;
  Index = g_new(gint,1);
  *Index = count;
  g_hash_table_insert(hash, head->indv, Index);
  head = head->next;  
  while(head != NULL)
    {
      if(g_hash_table_contains(hash,head->indv) == FALSE)
	{
	  count++;
	  Index = g_new(gint,1);
	  *Index = count;
	  g_hash_table_insert(hash, head->indv, Index);
	}
      else
	{
	  fprintf(stderr,"ERROR: The individual '%s' appears more than once\n",head->indv);
	  exit(1);
	}
      head = head->next;
    }
}
  
void missing_parent_check(family *head, GHashTable* hash)
{
  family *current;
  current = head;
  int count = 0;
  gint *Index;
  while(head != NULL)
    {
      count++;
      Index = g_new(gint,1);
      *Index = count;
      g_hash_table_insert(hash, head->indv, Index);
      head = head->next;
    }
  
  int lineNo = 0;

  while(current != NULL)
    {
      lineNo++;
      if (strcmp(current->dad,"0") != 0 && strcmp(current->mom,"0") != 0) /* Both parents NOT founder*/
      	{
  	  if(g_hash_table_contains(hash,current->dad) == FALSE || g_hash_table_contains(hash,current->mom) == FALSE)
  	    {
	      if (g_hash_table_contains(hash,current->dad) == FALSE)
		{
		  fprintf(stderr,"ERROR at Line %d: Missing information on Father %s\n",lineNo,current->dad);
		  exit(1);
		}
	      else if (g_hash_table_contains(hash,current->mom) == FALSE)
		{
		  fprintf(stderr,"ERROR at Line %d: Missing information on Mother %s\n",lineNo,current->mom);
		  exit(1);
		}
	      else if (g_hash_table_contains(hash,current->dad) == FALSE && g_hash_table_contains(hash,current->mom) == FALSE)
		{
		  fprintf(stderr,"ERROR at Line %d: Missing information on BOTH parents\n",lineNo);
		  exit(1);
		}
  	    }
  	}
      else if(strcmp(current->dad,"0") != 0 && strcmp(current->mom,"0") == 0) /* Dad is NOT a founder BUT Mom is a founder*/
	{
	  if (g_hash_table_contains(hash,current->dad) == FALSE)
	    {
	      fprintf(stderr,"ERROR at Line %d: Missing information on Father %s\n",lineNo,current->dad);
	      exit(1);
	    }
	}
      else if(strcmp(current->dad,"0") == 0 && strcmp(current->mom,"0") != 0) /* Dad is a founder BUT Mom is NOT a founder*/
	{
	  if (g_hash_table_contains(hash,current->mom) == FALSE)
	    {
	      fprintf(stderr,"ERROR at Line %d: Missing information on Mother %s\n",lineNo,current->mom);
	      exit(1);
	    }
	}
      current = current->next;
    }
}

void child_is_parent_of_parent(family *head) /* error check: Child is a parent of parent*/
{
  family *current;
  current = head;
  char C[MAX_NAME];
  char F[MAX_NAME];
  char M[MAX_NAME];
  int line_1 = 0;
  int line_2 = 0;
  
  while(head != NULL)
    {
      line_1++;
      line_2++;
      strcpy(C,head->indv);
      strcpy(F,head->dad);
      strcpy(M,head->mom);
      current = head->next;
      while (current != NULL)
	{
	  line_2++;
	  if((strcmp(F,current->indv) == 0 && strcmp(C,current->dad) == 0) || (strcmp(M,current->indv) == 0 && strcmp(C,current->mom) == 0))
	    {
	      fprintf(stderr,
		      "Data entry at Line %d: Child: %s\t Father: %s\t Mother: %s\n\n"
		      "ERROR at Line %d: Child is a parent of a parent (LOGICALLY IMPOSSIBLE)\n"
		      "Child: %s\t Father: %s\t Mother: %s\n",line_1,C,F,M,line_2,current->indv,current->dad,current->mom);
	      exit(1);
	    }
	  current = current->next;
	}
      head = head->next;
    }

}

void parent_of_ancestor_check(char *x, char *y, GHashTable *hash) /* error check: descendant is parent of an ancestor */
{
      if(never_visited == 0 && strcmp(x,y) == 0)
	{
	  fprintf(stderr,"ERROR: A descendant cannot be a parent of an ancestor (LOGICALLY IMPOSSIBLE)\n\n"
		  "Ancestor of a non-founder individual '%s' should always trace back to founders but NEVER itself\n",y);
	  exit(1);
	}

      else
	{
	  never_visited = 0;
	  if(strcmp(x,"0") != 0)
	    {
	      parent_of_ancestor_check(((guardian*)g_hash_table_lookup(hash,x))->DAD,y,hash);
	      parent_of_ancestor_check(((guardian*)g_hash_table_lookup(hash,x))->MOM,y,hash);
	    }
	}
      never_visited = 1;

}


void logical_error_check(family *head)
{
  GHashTable* hash_1 = g_hash_table_new_full(g_str_hash, g_str_equal, NULL, g_free); /* Creating Hash Table for missing parent error check */  
  missing_parent_check(head,hash_1);

  GHashTable* hash_2 = g_hash_table_new_full(g_str_hash, g_str_equal, NULL, g_free); /* Creating Hash Table for duplicate row error check */
  duplicate_row_check(head,hash_2);
  
  child_is_parent_of_parent(head); /* Special case: Child is a parent of a parent check */

  /* creating Hash Table for parent of ancestor check */
  family *current;
  current = head;
  guardian *P;
  GHashTable* hash_3 = g_hash_table_new_full(g_str_hash, g_str_equal,NULL, g_free); /* KEY: current->indv  ,  VALUE: Parent of current->indv */
  while(current != NULL)
    {
      P = g_new(guardian,1);
      strcpy(P->DAD,current->dad);
      strcpy(P->MOM,current->mom);
      g_hash_table_insert(hash_3,current->indv,P);
      current = current->next;
    }

  /* Running the parent of ancestor check over the linked list */
  current = head;
  while(current != NULL)
    {
      if(strcmp(((guardian*)g_hash_table_lookup(hash_3,current->indv))->DAD,"0") != 0 && strcmp(((guardian*)g_hash_table_lookup(hash_3,current->indv))->MOM,"0") != 0) /* We need to check the parent_of_ancestor only over the non-founder individuals */
	{
	  parent_of_ancestor_check(current->indv,current->indv,hash_3);
	}
      current = current->next;
    }

  g_hash_table_destroy(hash_1);
  g_hash_table_destroy(hash_2);
  g_hash_table_destroy(hash_3);  
}

void arrange(family *head, int **i2, int **f2, int **m2, int *n)
{
  family* current;

  GHashTable* hash = g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free); /* Creating Hash Table for all individuals */
  GHashTable* hash_founder = g_hash_table_new(g_str_hash, g_str_equal); /* Creating Hash Table for all founders */
  int number_of_founder,familyIndex;
  int *individual, *father, *mother; /* Creating pointer for the final ordered representation of the data */
  int count;

  count = node_count(head);

  Hash_Function_creation(head,hash,hash_founder, &number_of_founder, &familyIndex,count);
  printf("\n");
  individual = malloc(familyIndex * sizeof(int));
  father = malloc(familyIndex * sizeof(int));
  mother = malloc(familyIndex * sizeof(int));
  *i2 = malloc(familyIndex * sizeof(int));
  *f2 = malloc(familyIndex * sizeof(int));
  *m2 = malloc(familyIndex * sizeof(int));

  current = head;
  while(current != NULL)
    {
      int order;
      order = *(gint*)g_hash_table_lookup(hash,current->indv);
      if(order<=number_of_founder)
	{
	  *(individual + (order-1)) = order;
	  *(father + (order-1)) = 0;
	  *(mother + (order-1)) = 0;
	}
      else
	{
	  *(individual + (order-1)) = order;
	  *(father + (order-1)) = *(gint*)g_hash_table_lookup(hash,current->dad);
	  *(mother + (order-1)) = *(gint*)g_hash_table_lookup(hash,current->mom);
	}
      current = current->next;
    }
  for (int j = 0; j < familyIndex ; ++j)
    {
      *(*i2+j) = *(individual+j);
      *(*f2+j) = *(father+j);
      *(*m2+j) = *(mother+j);
    }
  current = head;
  while (current != NULL)
    {
      printf("Index: %d  Individual: %s\n",*(gint*)g_hash_table_lookup(hash,current->indv),current->indv);
      current = current->next;
    }
  printf("\n");
  g_hash_table_foreach(hash_founder, (GHFunc)iterator, "Index: %d is the Founder Individual: %s\n");
  printf("\n");
  printf("The indexed list of the pedigree arranged as follows :\n");
  for(int j = 0; j < familyIndex ; ++j)
    {
      printf("%d %d %d\n",*(*i2+j),*(*f2+j),*(*m2+j));
    }
  
  *n = familyIndex;
  g_hash_table_destroy(hash);
  g_hash_table_destroy(hash_founder);
  free(individual);
  free(father);
  free(mother);
}

