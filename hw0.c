#include <stdlib.h>
#include <stdio.h>
#include <time.h>

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
  float a = 0.05;
  float b = 0.1;
  float c = 0.4;
  float t = 0.1;
  //int n = 16384 + 2;
  int n = 5;

  int countX, countY = 0;

  // Allocate x and y to nxn
  float (*x)[n] = malloc(sizeof(*x) * n);
  float (*y)[n] = malloc(sizeof(*y) * n);

  initialize(n, x);
  smooth(n, x, y, a, b, c);
  count(n, x, t, &countX);
  count(n, y, t, &countY);

  //printArray(n, x);

  return 0;

}
