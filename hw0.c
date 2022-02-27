#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

void initialize(int n, float (*arr)[n]) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      arr[i][j] = random() / (float) RAND_MAX;
    }
  }
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
  for (int i = 1; i < n-1; i++) {
    for (int j = 1; j < n-1; j++) {
      B[i][j] = a*(A[i-1][j-1] + A[i-1][j+1] + A[i+1][j-1] + A[i+1][j+1]) +
        b*(A[i-1][j] + A[i+1][j] + A[i][j-1] + A[i][j+1]) +
        c*A[i][j];
    }
  }
}

void count(int n, float (*A)[n], float t, int* out) {
  int count = 0;

  for (int i = 1; i < n-1; i++) {
    for (int j = 1; j < n-1; j++) {
      if (A[i][j] < t) {
        count = count + 1;
      }
    }
  }
  *out = count;
}

int main() {
  /* Padding for printf */
  int pad1 = -40;
  int pad2 = 15;

  float a = 0.05;
  float b = 0.1;
  float c = 0.4;
  float t = 0.1;
  int n = 16384 + 2;
  //int n = 5;

  int countX, countY = 0;

  // Allocate x and y to nxn
  float (*x)[n] = malloc(sizeof(*x) * n);
  float (*y)[n] = malloc(sizeof(*y) * n);

  initialize(n, x);
  smooth(n, x, y, a, b, c);
  count(n, x, t, &countX);
  count(n, y, t, &countY);

  int elems = n*n;
  int inner_elems = (n-2)*(n-2);
  float fracX = countX / (double)inner_elems;
  float fracY = countY / (double)inner_elems;
  float arrSize = elems * 4.0 * 1e-9;

  /* Print info */
  printf("Summary\n");
  printf("-------\n");
  printf("%*s :: %*d\n", pad1, "Number of elements in a row/column", pad2, n);
  printf("%*s :: %*d\n", pad1, "Number of inner elements in a row/column", pad2, n-2);
  printf("%*s :: %*d\n", pad1, "Total number of elements", pad2, n*n);
  printf("%*s :: %*d\n", pad1, "Total number of inner elements", pad2, (n-2)*(n-2));
  printf("%*s :: %*.5f\n", pad1, "Memory (GB) used per array", pad2, arrSize);
  printf("%*s :: %*.2f\n", pad1, "Threshold", pad2, t);
  printf("%*s :: %.2f %.2f %.2f\n", pad1, "Smoothing constants (a,b,c)",a,b,c);
  printf("%*s :: %*d\n", pad1, "Number of elements below threshold (X)", pad2, countX);
  printf("%*s :: %*.5e\n", pad1, "Fraction of elements below threshold (X)", pad2, fracX);
  printf("%*s :: %*d\n", pad1, "Number of elements below threshold (Y)", pad2, countY);
  printf("%*s :: %*.5e\n", pad1, "Fraction of elements below threshold (Y)", pad2, fracY);

  return 0;

}
