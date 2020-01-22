#include <cstdio>
#include <cstdlib>
#include <cstddef>
#include <cstring>
#include <immintrin.h>
#include "ttmat.h"
#include <omp.h>

void loadTTMat(
    char *fileName,
    TTMat *mat)
{
  FILE *file = fopen(fileName, "rb");
  if (!file) { printf("ERROR: Unable to open the file %s.\n", fileName); exit(1); }
  fread(&(mat->d), sizeof(mat->d), 1, file);
  mat->m = (int *)malloc(mat->d * sizeof(mat->m[0])); // Number of matrix rows in each dimension
  fread(mat->m, sizeof(mat->m[0]), mat->d, file);
  mat->n = (int *)malloc(mat->d * sizeof(mat->n[0])); // Number of matrix columns in each dimension
  fread(mat->n, sizeof(mat->n[0]), mat->d, file);
  mat->r = (int *)malloc((mat->d + 1) * sizeof(mat->r[0]));
  mat->r[0] = mat->r[mat->d] = 1; // The first and the last ranks are always one.
  fread(mat->r + 1, sizeof(mat->r[0]), mat->d - 1, file); // Read the remaining D - 1 ranks
  mat->dimMatBegin = (size_t *)malloc((mat->d + 1) * sizeof(mat->dimMatBegin[0])); 
  mat->dimMatBegin[0] = 0;
  for (int d = 0; d < mat->d; d++) { // Setup pointers for the beginning of data for the matrix of each dimension
    // Each dimension has a matrix of size m[d] x n[d], each element of this matrix is of size r[d] x r[d + 1].
    mat->dimMatBegin[d + 1] = mat->dimMatBegin[d] + mat->m[d] * mat->n[d] * mat->r[d] * mat->r[d + 1];
  }
  mat->data = (double *)malloc(mat->dimMatBegin[mat->d] * sizeof(mat->data[0])); // Allocate all matrices
  fread(mat->data, sizeof(mat->data[0]), mat->dimMatBegin[mat->d], file); // Read the input
  fclose(file);
}

void destroyTTMat(TTMat *mat)
{
  free(mat->m);
  free(mat->n);
  free(mat->r);
  free(mat->dimMatBegin);
  free(mat->data);
}

void printTTMat(
    TTMat *mat,
    FILE *outFile)
{
  fprintf(outFile, "TTMat:\n");
  fprintf(outFile, "D = %d\n", mat->d);
  for (int d = 0; d < mat->d; d++) {
    fprintf(outFile, "M[%d]xN[%d]=%dx%d\t\t", d, d, mat->m[d], mat->n[d]);
    fprintf(outFile, "R[%d]xR[%d]=%dx%d\n", d, d + 1, mat->r[d], mat->r[d + 1]);
  }
  for (int d = 0; d < mat->d; d++) {
    fprintf(outFile, "Printing matrices for dimension %d:\n", d);
    for (int i = 0; i < mat->m[d]; i++) {
      for (int j = 0; j < mat->n[d]; j++) {
        fprintf(outFile, "  Printing the matrix block (%d,%d):\n", i, j);
        double *blockBegin = mat->data + mat->dimMatBegin[d] + (i + j * mat->m[d]) * (mat->r[d] * mat->r[d + 1]);
        for (int ri = 0; ri < mat->r[d]; ri++) {
          for (int rj = 0; rj < mat->r[d + 1]; rj++) {
            fprintf(outFile, "    %lf ", blockBegin[ri + rj * mat->r[d]]);
          }
          fprintf(outFile, "\n");
        }
      }
    }
  }
}

double *getTTMatBlock(
    TTMat *mat,
    int dim,
    int rowIdx,
    int colIdx)
{
  return mat->data + mat->dimMatBegin[dim] + ((size_t)(rowIdx + colIdx * mat->m[dim]) * (size_t)(mat->r[dim] *
        mat->r[dim + 1]));
}

void multiplyAddKronecker(
    double *a,
    int anrows,
    int ancols,
    double *b,
    int bnrows,
    int bncols,
    double *c)
{
  size_t  ka = (size_t)(anrows * bnrows);
  size_t kb = (size_t)(bnrows * bncols);

  for (int ja = 0; ja < ancols; ja++) {
    int a_remaining = anrows%4;
    for (int ia = 0; ia < anrows-a_remaining; ia++) {
      size_t na = (ia + ja * anrows);
      size_t la = na*kb;
      double aelem = a[na];
      
      for (int jb = 0; jb < bncols; jb++) {
	size_t  lb = jb * bnrows;
	size_t ma = jb *ka;

	int b_remaining = bnrows%4;
	for (int ib = 0; ib < bnrows-b_remaining; ib++) {
          double belem = b[ib + lb];
          size_t cElemIdx = la + ib + ma;
          c[cElemIdx] += aelem * belem;
	}

        for (int ib = bnrows-b_remaining; ib < bnrows; ib++) {
          double belem = b[ib + lb];
          size_t cElemIdx = la + ib + ma;
          c[cElemIdx] += aelem * belem;
        }
      }
    }
    for (int ia = anrows-a_remaining; ia < anrows; ia++) {
      size_t na = (ia + ja * anrows);
      size_t la = na*kb;
      double aelem = a[na];
      for (int jb = 0; jb < bncols; jb++) {
	size_t  lb = jb * bnrows;
	size_t ma = jb *ka;
        for (int ib = 0; ib < bnrows; ib++) {
          double belem = b[ib + lb];
          size_t cElemIdx = la + ib + ma;
          c[cElemIdx] += aelem * belem;
        }
      }
    }
  }
}

void multiplyTTMatVec(
    TTMat *A,
    TTVec *x,
    TTVec *y)
{
  // Allocate y
  y->d = x->d;
  y->m = (int *)malloc(y->d * sizeof(y->m[0]));
  y->r = (int *)malloc((y->d + 1) * sizeof(y->r[0]));
  y->dimVecBegin = (size_t *)malloc((y->d + 1) * sizeof(y->dimVecBegin[0]));
  y->dimVecBegin[0] = 0;
  y->r[0] = A->r[0] * x->r[0];
  for (int d = 0; d < y->d; d++) {
    y->m[d] = A->m[d];
    y->r[d + 1] = A->r[d + 1] * x->r[d + 1];
    y->dimVecBegin[d + 1] = y->dimVecBegin[d] + y->m[d] * y->r[d] * y->r[d + 1];
  }
  y->data = (double *)malloc(y->dimVecBegin[y->d] * sizeof(y->data[0]));
  memset(y->data, 0, y->dimVecBegin[y->d] * sizeof(y->data[0]));

  // Now perform the matrix-vector multiplication in each dimension
  for (int d = 0; d < y->d; d++) {
    // for each matrix at line d of A 
    for (int m = 0; m < A->m[d]; m++) {
      double *ymBlockBegin = getTTVecBlock(y, d, m);
      // for each matrix at column d of A
      for (int n = 0; n < A->n[d]; n++) {
        double *AmnBlockBegin = getTTMatBlock(A, d, m, n);
        double *xnBlockBegin = getTTVecBlock(x, d, n);
        multiplyAddKronecker(AmnBlockBegin, A->r[d], A->r[d + 1], xnBlockBegin, x->r[d], x->r[d + 1], ymBlockBegin);
      }
    }
  }
}
 

//y->d * A->m[d] * A->n[d] * ancols * anrows * bncols * bnrows * 2
//D    * N       * N       * Ra     * Qa     * Rx     * Qx     * 2
