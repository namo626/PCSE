#include <stdlib.h>
#include <stdio.h>

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
}

int main() {
  float a = 0.05;
  float b = 0.1;
  float c = 0.4;
  float t = 0.1;
  //int n = 16384 + 2;
  int n = 5;

  // Allocate x and y to nxn
  float (*x)[n] = malloc(sizeof(*x) * n);
  float (*y)[n] = malloc(sizeof(*y) * n);

  initialize(n, x);
  initialize(n, y);

  printArray(n, x);

  return 0;

}
