#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

void initialize(int rank, int NP, int m, int n, float (*arr)[n]) {
  if (rank == 0) {
    for (int k = 1; k <= NP; k++) {
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          arr[i][j] = random() / (float) RAND_MAX;
        }
      }
      if (k < NP) {
        MPI_Send(&arr[0][0], m*n, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
      }
    }
  } else {
    MPI_Recv(&arr[0][0], m*n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

/* Count CPU time in sec */
float getTime(clock_t i1, clock_t i2) {
  float t1 = (i2 - i1) / (float)CLOCKS_PER_SEC;
  return t1;
}

/* Get the neighbors */
void getNeighbors(int rank, int M, int N, int* north, int* south, int* east, int* west) {
  int total = M * N;
  int n0, s0, e0, w0; // dummy vars

  n0 = rank + N;
  // if neighbor exceeds limit, this tile is at a boundary
  if (n0 > total - 1) {
    *north = -1;
  } else {
    *north = n0;
  }

  s0 = rank - N;
  if (s0 < 0) {
    *south = -1;
  } else {
    *south = s0;
  }

  e0 = rank + 1;
  if (e0 % N != 0) {
    *east = e0;
  } else {
    *east = -1;
  }

  w0 = rank - 1;
  if (rank % N != 0) {
    *west = w0;
  } else {
    *west = -1;
  }
}

int main() {
  MPI_Init(NULL, NULL);

  /* Size of local array mxn */
  int m = 500;
  int n = 400;

  /* Get the number of processes (same as number of tiles), MxN
     Assume that M x N = NP
  */
  int NP;
  MPI_Comm_size(MPI_COMM_WORLD, &NP);

  int M = 2;
  int N = 2;

  // Get the rank of the process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int north, south, east, west;
  getNeighbors(rank, M, N, &north, &south, &east, &west);
  //printf("PE: %d, Neighbors = %d, %d, %d, %d\n", rank, west, north, east, south);

  /* Initialize local array mxn */
  float (*x)[n] = malloc(sizeof(*x) * m);
  initialize(rank, NP, m, n, x);
  printf("PE: %d, x[m-1][n-1] = %.4f\n", rank, x[m-1][n-1]);

  /* Assign rank to tile */


  /* /\* Timing variables *\/ */
  /* clock_t i1, i2; */
  /* float t1, t2, t3, t4, t5, t6; */

  /* float a = 0.05; */
  /* float b = 0.1; */
  /* float c = 0.4; */
  /* float t = 0.1; */
  /* int n = 16384 + 2; */
  /* //int n = 5; */

  /* int countX, countY = 0; */

  /* // Allocate x and y to nxn */
  /* i1 = clock(); */
  /* float (*x)[n] = malloc(sizeof(*x) * n); */
  /* i2 = clock(); */
  /* t1 = getTime(i1, i2); */

  /* i1 = clock(); */
  /* float (*y)[n] = malloc(sizeof(*y) * n); */
  /* i2 = clock(); */
  /* t2 = getTime(i1, i2); */

  /* i1 = clock(); */
  /* initialize(n, x); */
  /* i2 = clock(); */
  /* t3 = getTime(i1, i2); */

  /* i1 = clock(); */
  /* smooth(n, x, y, a, b, c); */
  /* i2 = clock(); */
  /* t4 = getTime(i1, i2); */

  /* i1 = clock(); */
  /* count(n, x, t, &countX); */
  /* i2 = clock(); */
  /* t5 = getTime(i1, i2); */

  /* i1 = clock(); */
  /* count(n, y, t, &countY); */
  /* i2 = clock(); */
  /* t6 = getTime(i1, i2); */

  /* int elems = n*n; */
  /* int inner_elems = (n-2)*(n-2); */
  /* float fracX = countX / (double)inner_elems; */
  /* float fracY = countY / (double)inner_elems; */
  /* float arrSize = elems * 4.0 * 1e-9; */

  /* Padding for printf */
  int pad1 = -40;
  int pad2 = 15;
  int pad3 = -20;
  int pad4 = 10;

  /* /\* Print info *\/ */
  /* printf("Summary\n"); */
  /* printf("-------\n"); */
  /* printf("%*s :: %*d\n", pad1, "Number of elements in a row/column", pad2, n); */
  /* printf("%*s :: %*d\n", pad1, "Number of inner elements in a row/column", pad2, n-2); */
  /* printf("%*s :: %*d\n", pad1, "Total number of elements", pad2, n*n); */
  /* printf("%*s :: %*d\n", pad1, "Total number of inner elements", pad2, (n-2)*(n-2)); */
  /* printf("%*s :: %*.5f\n", pad1, "Memory (GB) used per array", pad2, arrSize); */
  /* printf("%*s :: %*.2f\n", pad1, "Threshold", pad2, t); */
  /* printf("%*s :: %.2f %.2f %.2f\n", pad1, "Smoothing constants (a,b,c)",a,b,c); */
  /* printf("%*s :: %*d\n", pad1, "Number of elements below threshold (X)", pad2, countX); */
  /* printf("%*s :: %*.5e\n", pad1, "Fraction of elements below threshold (X)", pad2, fracX); */
  /* printf("%*s :: %*d\n", pad1, "Number of elements below threshold (Y)", pad2, countY); */
  /* printf("%*s :: %*.5e\n", pad1, "Fraction of elements below threshold (Y)", pad2, fracY); */

  /* printf("\n"); */

  /* printf("Action\n"); */
  /* printf("------\n"); */
  /* printf("%*s :: %*.3f\n", pad3, "CPU: Alloc-X", pad4, t1); */
  /* printf("%*s :: %*.3f\n", pad3, "CPU: Alloc-Y", pad4, t2); */
  /* printf("%*s :: %*.3f\n", pad3, "CPU: Init-X", pad4, t3); */
  /* printf("%*s :: %*.3f\n", pad3, "CPU: Smooth", pad4, t4); */
  /* printf("%*s :: %*.3f\n", pad3, "CPU: Count-X", pad4, t5); */
  /* printf("%*s :: %*.3f\n", pad3, "CPU: Count-Y", pad4, t6); */

  MPI_Finalize();
  return 0;

}
