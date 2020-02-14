#include <uGnix.h>

typedef struct ancinter
{
  double lower;
  double upper;
} anc_interval;



anc_interval* intersect(anc_interval* ancInt1, anc_interval* ancInt2)
{
  double newLower = MAX(ancInt1->lower,ancInt2->lower);
  double newUpper = MIN(ancInt1->upper,ancInt2->upper);
  anc_interval* intersection;
  if(newUpper <= newLower)
    return NULL;
  else
    {
      intersection = malloc(sizeof(anc_interval));
      intersection->lower = newLower;
      intersection->upper = newUpper;
      return intersection;
    }
}

char version[] = "coalsim";

int main()
{
  anc_interval* inter1;
  anc_interval* inter2;
  anc_interval* intersect12;

  inter1 = malloc(sizeof(anc_interval));
  inter2 = malloc(sizeof(anc_interval));
  intersect12 = malloc(sizeof(anc_interval));
  inter1->lower = 0.3;
  inter1->upper = 0.6;
  inter2->lower = 0.2;
  inter2->upper = 0.9;
  intersect12 = intersect(inter1,inter2);
  if(intersect12!=NULL)
    printf("intersection is: (%f, %f)\n",intersect12->lower,intersect12->upper);
  else
    printf("intersection is empty set\n");
}

