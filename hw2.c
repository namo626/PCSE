#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>


/* m and n include the ghost cells, so don't update those */
void initialize(int rank, int NP, int m, int n, float (*arr)[n]) {
  if (rank == 0) {
    for (int k = 1; k <= NP; k++) {
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
          if (i == 0 || i == m-1 || j == 0 || j == n-1) {
            arr[i][j] = 0.0;
          } else {
            arr[i][j] = random() / (float) RAND_MAX;
          }
          //arr[i][j] = 1.0;
        }
      }
      if (k < NP) {
        MPI_Send(&arr[0][0], m*n, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
      }
    }
  } else {
    MPI_Recv(&arr[0][0], m*n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  MPI_Barrier(MPI_COMM_WORLD);

}

void printArray(int m, int n, float (*arr)[n]) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("%.3f ", arr[i][j]);
    }
    printf("\n");
  }
}

void setZero(int mm, int nn, float (*arr)[nn]) {
  for (int i = 0; i < mm; i++) {
    for (int j = 0; j < nn; j++) {
      arr[i][j] = 0.0;
    }
  }
}

void smooth(int mm, int nn, float (*A)[nn], float (*B)[nn], float b, float c, int north, int south, int east, int west) {

  int m = mm - 2;
  int n = nn - 2;

  /* Determine the range based on neighbor */
  int imin = (north != -1) ? 1 : 2;
  int imax = (south != -1) ? m : m-1;
  int jmin = (west != -1) ? 1 : 2;
  int jmax = (east != -1) ? n : n-1;

  for (int i = imin; i <= imax; i++) {
    for (int j = jmin; j <= jmax; j++) {
      B[i][j] = b*(A[i-1][j] + A[i+1][j] + A[i][j-1] + A[i][j+1]) +
        c*A[i][j];
    }
  }
}

void count(int mm, int nn, float (*arr)[nn], float t, int north, int south, int east, int west, int* out) {
  int count = 0;
  int m = mm - 2;
  int n = nn - 2;

  /* Determine the range based on neighbor */
  int imin = (north != -1) ? 1 : 2;
  int imax = (south != -1) ? m : m-1;
  int jmin = (west != -1) ? 1 : 2;
  int jmax = (east != -1) ? n : n-1;

  for (int i = imin; i <= imax; i++) {
    for (int j = jmin; j <= jmax; j++) {
      if (arr[i][j] < t) {
        count = count + 1;
      }
    }
  }
  *out = count;
}

/* Copy the values at ghost cells */
void updateCells(MPI_Datatype Column, int mm, int nn, float (*arr)[nn], int north, int south, int east, int west) {
  MPI_Request reqs[4];
  MPI_Status stats[4];
  MPI_Request dummy[4];

  int m = mm - 2;
  int n = nn - 2;

  if (north != -1) {
    /* copy the bottom ghost row of the north tile */
    MPI_Isend(&arr[1][1], n, MPI_FLOAT, north, 0, MPI_COMM_WORLD, &dummy[0]);
    MPI_Irecv(&arr[0][1], n, MPI_FLOAT, north, 0, MPI_COMM_WORLD, &reqs[0]);
  } else {
    reqs[0] = MPI_REQUEST_NULL;
  }

  if (south != -1) {
    MPI_Isend(&arr[m][1], n, MPI_FLOAT, south, 0, MPI_COMM_WORLD, &dummy[1]);
    MPI_Irecv(&arr[m+1][1], n, MPI_FLOAT, south, 0, MPI_COMM_WORLD, &reqs[1]);
  } else {
    reqs[1] = MPI_REQUEST_NULL;
  }

  if (east != -1) {
    MPI_Isend(&arr[1][n], 1, Column, east, 0, MPI_COMM_WORLD, &dummy[2]);
    MPI_Irecv(&arr[1][n+1], 1, Column, east, 0, MPI_COMM_WORLD, &reqs[2]);
  } else {
    reqs[2] = MPI_REQUEST_NULL;
  }

  if (west != -1) {
    MPI_Isend(&arr[1][1], 1, Column, west, 0, MPI_COMM_WORLD, &dummy[3]);
    MPI_Irecv(&arr[1][0], 1, Column, west, 0, MPI_COMM_WORLD, &reqs[3]);
  } else {
    reqs[3] = MPI_REQUEST_NULL;
  }

  // Wait for all ghost cells to be updated
  MPI_Waitall(4, reqs, stats);
}

/* Compute the closest pair of factors of a given number P */
void closestFactors(int P, int* M, int* N) {
  int c = (int) sqrt(P);
  while (P % c != 0) {
    c = c - 1;
  }
  *M = c;
  *N = P / c;
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
  *north = (n0 > total-1) ? -1 : n0;

  s0 = rank - N;
  *south = (s0 < 0) ? -1 : s0;

  e0 = rank + 1;
  *east = (e0 % N != 0) ? e0 : -1;

  w0 = rank - 1;
  *west = (rank % N != 0) ? w0 : -1;
}

int main() {
  MPI_Init(NULL, NULL);

  /* Size of local array mxn, fixed for every processor */
  int m = 1500;
  int n = 1600;

  int NP;
  MPI_Comm_size(MPI_COMM_WORLD, &NP);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Number of rows and number of columns; NP = M * N */
  int M, N;
  closestFactors(NP, &M, &N);
  if (rank == 0) {
    printf("M = %d, N = %d\n", M, N);
  }

  /* Compute the neighbors */
  int north, south, east, west;
  getNeighbors(rank, M, N, &north, &south, &east, &west);
  printf("PE: %d, Neighbors = %d, %d, %d, %d\n", rank, west, north, east, south);

  /* Timing variables */
  clock_t i1, i2;
  float t1, t2, t3, t4, t5, t6;

  i1 = clock();
  /* Initialize local array mxn. Allocate extra space for ghost cells */
  float (*x)[n+2] = malloc(sizeof(*x) * (m+2));
  i2 = clock();
  t1 = getTime(i1, i2);

  i1 = clock();
  float (*y)[n+2] = malloc(sizeof(*y) * (m+2));
  i2 = clock();
  t2 = getTime(i1, i2);
  setZero(m+2, n+2, y);

  /* Initialize includes sharing to other processors */
  i1 = clock();
  initialize(rank, NP, m+2, n+2, x);
  i2 = clock();
  t3 = getTime(i1, i2);

  /* Data type for columns */
  MPI_Datatype Column;
  MPI_Type_vector(m, 1, n+2, MPI_FLOAT, &Column);
  MPI_Type_commit(&Column);

  /* if (rank == 0) { */
  /*   printf("Rank 0 before update:\n"); */
  /*   printArray(m+2, n+2, x); */
  /*   printf("\n"); */
  /* } */
  /* Update ghost cells */
  updateCells(Column, m+2, n+2, x, north, south, east, west);

  /* Smoothing */
  float b = 0.1;
  float c = 0.4;
  float t = 0.1;

  i1 = clock();
  smooth(m+2, n+2, x, y, b, c, north, south, east, west);
  i2 = clock();
  t4 = getTime(i1, i2);

  int counts[2] = {0, 0};
  int globalCounts[2] = {0, 0};

  i1 = clock();
  count(m+2, n+2, x, t, north, south, east, west, &counts[0]);
  i2 = clock();
  t5 = getTime(i1, i2);

  i1 = clock();
  count(m+2, n+2, y, t, north, south, east, west, &counts[1]);
  i2 = clock();
  t6 = getTime(i1, i2);

  MPI_Reduce(counts, globalCounts, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  // testing
  float (*z)[n+2] = malloc(sizeof(*z) * (m+2));

  /* if (rank == 0) { */
  /*   MPI_Recv(&z[0][0], (m+2)*(n+2), MPI_FLOAT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); */
  /*   printf("Rank 1 X:\n"); */
  /*   printArray(m+2, n+2, z); */
  /*   printf("\n"); */


  /*   printf("Rank 0 X:\n"); */
  /*   printArray(m+2, n+2, x); */
  /*   printf("\n"); */

  /*   printf("Rank 0 Y:\n"); */
  /*   printArray(m+2, n+2, y); */
  /*   printf("\n"); */
  /* } */

  /* if (rank == 1) { */
  /*   MPI_Send(&x[0][0], (m+2)*(n+2), MPI_FLOAT, 0, 0, MPI_COMM_WORLD); */
  /* } */


  int globalElems = M*m*N*n;
  int inner_elems = (M*m-2)*(N*n-2);
  float fracX = globalCounts[0] / (double)inner_elems;
  float fracY = globalCounts[1] / (double)inner_elems;
  float arrSize = (m+2)*(n+2) * 4.0 * 1e-9;

  /* Padding for printf */
  int pad1 = -40;
  int pad2 = 15;
  int pad3 = -20;
  int pad4 = 10;

  if (rank == 0) {
  /* Print info */
    printf("Summary\n");
    printf("-------\n");
  /* printf("%*s :: %*d\n", pad1, "Number of elements in a row/column", pad2, n); */
  /* printf("%*s :: %*d\n", pad1, "Number of inner elements in a row/column", pad2, n-2); */
  /* printf("%*s :: %*d\n", pad1, "Total number of elements", pad2, n*n); */
  /* printf("%*s :: %*d\n", pad1, "Total number of inner elements", pad2, (n-2)*(n-2)); */
    printf("%*s :: %*.5f\n", pad1, "Memory (GB) used per local array", pad2, arrSize);
  /* printf("%*s :: %*.2f\n", pad1, "Threshold", pad2, t); */
  /* printf("%*s :: %.2f %.2f %.2f\n", pad1, "Smoothing constants (a,b,c)",a,b,c); */
    printf("%*s :: %*d\n", pad1, "Number of elements below threshold (X)", pad2, globalCounts[0]);
    printf("%*s :: %*.5e\n", pad1, "Fraction of elements below threshold (X)", pad2, fracX);
    printf("%*s :: %*d\n", pad1, "Number of elements below threshold (Y)", pad2, globalCounts[1]);
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
  }

  MPI_Finalize();
  return 0;

}
