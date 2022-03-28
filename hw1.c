#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "omp.h"

void initialize(int n, float (*arr)[n]) {
#ifdef _OPENMP
  #pragma omp parallel
  {
    unsigned int seed = omp_get_thread_num();
    #pragma omp for
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	arr[i][j] = rand_r(&seed) / (float) RAND_MAX;
      }
    }
  }
#else
  int seed = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      //arr[i][j] = random() / (float) RAND_MAX;
      arr[i][j] = rand_r(&seed) / (float) RAND_MAX;
    }
  }
#endif
}

void printArray(int n, float (*arr)[n]) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      printf("%.2f ", arr[i][j]);
    }
    printf("\n");
  }
}

void smooth(int n, float (*A)[n], float (*B)[n], float a, float b, float c) {
  #pragma omp parallel for schedule(static)
  for (int i = 1; i < n-1; i++) {
    for (int j = 1; j < n-1; j++) {
      B[i][j] = a*(A[i-1][j-1] + A[i-1][j+1] + A[i+1][j-1] + A[i+1][j+1]) +
        b*(A[i-1][j] + A[i+1][j] + A[i][j-1] + A[i][j+1]) +
        c*A[i][j];
    }
  }
}

void count(int n, float (*A)[n], float t, long* out) {
  long count = 0;

#pragma omp parallel for reduction(+:count)
  for (int i = 1; i < n-1; i++) {
    for (int j = 1; j < n-1; j++) {
      if (A[i][j] < t) {
        count = count + 1;
      }
    }
  }
  *out = count;
}

/* Choosing between serial or parallel */
double timeit() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return (double) clock() / CLOCKS_PER_SEC ;
#endif
}

/* Count CPU time in sec */
double getTime(double i1, double i2) {
  double t1 = (i2 - i1);
  return t1;
}

int main() {
  /* Padding for printf */
  int pad1 = -40;
  int pad2 = 20;
  int pad3 = -20;
  int pad4 = 10;

  /* Timing variables */
  double i1, i2;
  double t1, t2, t3, t4, t5, t6;

  float a = 0.05;
  float b = 0.1;
  float c = 0.4;
  float t = 0.1;
  //int n = 16384 + 2;
  long n = 98304+2;

  long countX, countY = 0;

  int num_threads = 1;

  /* Dummy parallel region */
#ifdef _OPENMP
  #pragma omp parallel
  {
  printf("This is thread %d\n", omp_get_thread_num());
#pragma omp single
  num_threads = omp_get_num_threads();
  }
#endif

  // Allocate x and y to nxn
  i1 = timeit();
  float (*x)[n] = malloc(sizeof(*x) * n);
  i2 = timeit();
  t1 = getTime(i1, i2);

  i1 = timeit();
  float (*y)[n] = malloc(sizeof(*y) * n);
  i2 = timeit();
  t2 = getTime(i1, i2);

  i1 = timeit();
  initialize(n, x);
  i2 = timeit();
  t3 = getTime(i1, i2);

  /* Problem 2b, repeat initialization */
  // printf("Init 1 : %.3fs\n", t3);

  // i1 = timeit();
  // initialize(n, x);
  // i2 = timeit();
  // t3 = getTime(i1, i2);

  // printf("Init 2 : %.3fs\n", t3);

  /* End problem 2b */

  i1 = timeit();
  smooth(n, x, y, a, b, c);
  i2 = timeit();
  t4 = getTime(i1, i2);

  i1 = timeit();
  count(n, x, t, &countX);
  i2 = timeit();
  t5 = getTime(i1, i2);

  i1 = timeit();
  count(n, y, t, &countY);
  i2 = timeit();
  t6 = getTime(i1, i2);

  long elems = n*n;
  long inner_elems = (n-2)*(n-2);
  float fracX = countX / (double)inner_elems;
  float fracY = countY / (double)inner_elems;
  float arrSize = elems * 4.0 * 1e-9;

  /* Print info */
  printf("Summary (Number of threads = %d)\n", num_threads);
  printf("-------\n");
  printf("%*s :: %*d\n", pad1, "Number of elements in a row/column", pad2, n);
  printf("%*s :: %*d\n", pad1, "Number of inner elements in a row/column", pad2, n-2);
  printf("%*s :: %*ld\n", pad1, "Total number of elements", pad2, elems);
  printf("%*s :: %*ld\n", pad1, "Total number of inner elements", pad2, inner_elems);
  printf("%*s :: %*.5f\n", pad1, "Memory (GB) used per array", pad2, arrSize);
  printf("%*s :: %*.2f\n", pad1, "Threshold", pad2, t);
  printf("%*s :: %.2f %.2f %.2f\n", pad1, "Smoothing constants (a,b,c)",a,b,c);
  printf("%*s :: %*d\n", pad1, "Number of elements below threshold (X)", pad2, countX);
  printf("%*s :: %*.5e\n", pad1, "Fraction of elements below threshold (X)", pad2, fracX);
  printf("%*s :: %*d\n", pad1, "Number of elements below threshold (Y)", pad2, countY);
  printf("%*s :: %*.5e\n", pad1, "Fraction of elements below threshold (Y)", pad2, fracY);

  printf("\n");

  printf("Action\n");
  printf("------\n");
  printf("%*s :: %*.3f\n", pad3, "CPU: Alloc-X", pad4, t1);
  printf("%*s :: %*.3f\n", pad3, "CPU: Alloc-Y", pad4, t2);
  printf("%*s :: %*.3f\n", pad3, "CPU: Init-X", pad4, t3);
  printf("%*s :: %*.3f\n", pad3, "CPU: Smooth", pad4, t4);
  printf("%*s :: %*.3f\n", pad3, "CPU: Count-X", pad4, t5);
  printf("%*s :: %*.3f\n", pad3, "CPU: Count-Y", pad4, t6);

  return 0;

}
