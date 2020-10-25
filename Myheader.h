#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
///llibreries que s'utilitzen



///declaració de totes les funcions que s'utilitzen en el programa
void previo(double *L, double *dx, double *a, int *N, int *k2,int *k1);
void introduce(int N,int n,int k1,int k2,int A,double **m1);
void triangular(int N,int n,double **m1);
void solution(int N,int n,double dx,double **m1, double *sol);


