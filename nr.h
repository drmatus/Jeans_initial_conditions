#include <math.h>
#include "nrutil.h"

#define TINY 1.0e-25 //A small number.
#define FREERETURN {free_vector(d,1,n);free_vector(c,1,n);return;}

void polint (float xa[], float ya[], int n, float x, float *y, float *dy);

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);

void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

void ratint(float xa[], float ya[], int n, float x, float *y, float *dy);
