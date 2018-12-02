#ifndef TTVEC_H
#define TTVEC_H

typedef struct
{
  int d;
  int *m, *r;
  double *data;
  size_t *dimVecBegin;
} TTVec;

void loadTTVec(
    char *fileName,
    TTVec *vec);

void saveTTVec(
    char *fileName,
    TTVec *vec);

void printTTVec(
    TTVec *vec,
    FILE *outFile);

double *getTTVecBlock(
    TTVec *vec,
    int dim,
    int rowIdx);

#endif