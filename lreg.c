/* @author: gregorydhill */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float** matrix(int n, int m) {
  float **X;
  X = malloc(n * sizeof (*X));
  for (int i=0; i<n; i++)
    X[i] = malloc(m * sizeof(*X[i]));
  return X;
}

void clear(int n, float** X){
  for (int i=0; i<n; i++)
    free(X[i]);
  free(X);
}

float** transpose(int n, int m, float** X) {
  float** X_ = matrix(m,n);
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      X_[i][j] = X[j][i];
    }
  }
  return X_;
}

float** product(int n, int m, int p, int q, float** A, float** B) {
  float** C = matrix(n,q);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < q; j++)
      C[i][j] = 0;

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < q; j++) {
      for (size_t k = 0; k < p; k++) {
        C[i][j] += A[i][k]*B[k][j];
      }
    }
  }

  return C;
}

void identity(int n, int m, float** X) {
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      if (i==j) X[i][j]=1;
      else X[i][j]=0;
    }
  }
}

float** get_minor(int row, int col, int n, float** M) {
  int k = 0;
  int l = 0;
  int s = n-1;
  float** m = matrix(s,s);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i!=row&&j!=col) {
        m[k][l] = M[i][j]; l++;
      }
    }
    if (i!=row) k++;
    l=0;
  }
  return m;
}

float determinant(int n, float** M) {
  if (n==2)
    return (M[0][0]*M[1][1])-(M[0][1]*M[1][0]);

  float det = 0;
  for (size_t j = 0; j < n; j++) {
    float** m = get_minor(0, j, n, M);
    det +=  pow((-1),j)*M[0][j]*determinant(n-1, m);
    clear(n, m);
  }
  return det;
}

float** inverse(int n, float** M) {
  float** C = matrix(n, n);
  float d = determinant(n, M);

  if (n==2) {
    C[0][0] = M[1][1]/d;
    C[0][1] = (-1)*M[0][1]/d;
    C[1][0] = (-1)*M[1][0]/d;
    C[1][1] = M[0][0]/d;
    return C;
  }

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      float** m = get_minor(i, j, n, M);
      C[i][j] = (pow((-1),i+j)*determinant(n-1, m))/d;
    }
  }

  float** A = transpose(n, n, C);
  clear(n, C);
  return A;
}

void output(int n, int m, float** X, char* T) {
  printf("%s\n---\n", T);
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < m; j++) {
      if (j==m-1) printf("%.2f", X[i][j]);
      else printf("%.2f, ", X[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

/*
* Linear Least Squares Optimization - Normal Equations Approach
*
* @param n: number of input rows
* @param m: number of input columns
* @param X: input matrix
* @param y: output vector
*
* @return: linear weights
*/
float** lstsq(int n, int m, float** X, float** y) {

  // TODO: include bias column

  float** X_ = transpose(n, m, X);
  float** A = product(m, n, n, m, X_, X);

  printf("\n");
  output(n, 1, y, "y");
  output(n, m, X, "X");
  output(m, n, X_, "X^T");
  output(m, m, A, "A=X*X^T");

  // non-square matrices do not have determinant
  float d = determinant(m, A);
  if (d==0) {
    printf("Matrix non-invertible.\n\n");
    exit(-1);
  }
  else printf("Determinant: %.1f\n\n", d);

  float** A_ = inverse(m, A);
  float** B = product(m, m, m, n, A_, X_);
  float** w = product(m, n, n, 1, B, y);

  output(m, m, A_, "A^T");
  output(m, n, B, "B=A^T*X^T");
  output(m, 1, w, "B*y");

  clear(n, X_);
  clear(n, A);
  clear(n, A_);
  clear(n, B);

  return w;
}

int main() {

  int n = 3;  // rows
  int m = 2;  // columns

  float** X = matrix(n,m);  // inputs
  float** y = matrix(n,1);  // outputs

  // test input
  X[0][0] = 1;
  X[0][1] = 3;
  X[1][0] = 2;
  X[1][1] = 4;
  X[2][0] = 1;
  X[2][1] = 6;

  // test output
  y[0][0] = 4;
  y[1][0] = 1;
  y[2][0] = 3;

  float** w = lstsq(n, m, X, y);  // weights

  float** f = product(n, m, m, 1, X, w);  // function values
  output(n, 1, f, "F=X*w");

  clear(n, X);
  clear(n, y);
  clear(n, w);
  clear(n, f);

  return 0;
}
