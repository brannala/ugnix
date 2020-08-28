#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<glib.h>
#include<regex.h>
#include<limits.h>
#include<unistd.h>
#include<ctype.h>
#include<stdbool.h>

#ifndef DATA_H_
#define DATA_H_

#define MAX_FILE_NAME 100
#define MAX_NAME 30

typedef struct trio
{
  char indv[MAX_NAME];
  char dad[MAX_NAME];
  char mom[MAX_NAME];
  struct trio* next;
} family;

typedef struct P
{
  char DAD[MAX_NAME];
  char MOM[MAX_NAME];
} guardian;


void  format_error_check(FILE *f, int *e); /* Checking data is in correct format for every trio */
void arrange(family *head, int **i2, int **f2, int **m2, int *n); /* Arrange data in the format required for Kinship Calculation*/
family* makeFamily(char* indv, char* dad, char* mom);
family* Linked_List_Creation(FILE *f);
void logical_error_check(family *head);
void parent_of_ancestor_check(char *x, char *y, GHashTable *hash);
void child_is_parent_of_parent(family *head);
void missing_parent_check(family *head,GHashTable* hash);
void duplicate_row_check(family *head, GHashTable *hash);
int node_count(family *head);
void Hash_Function_creation(family *head, GHashTable* hash, GHashTable* hash_founder, int *n_founder, int *n_family_index, int node_count );

#endif
