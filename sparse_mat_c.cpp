#include <cstdlib>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;
#include "Vector.hpp"
#include "full_mat_c.hpp"
#include "abstract_mat_c.hpp"
#include "sparse_mat_c.hpp"
#include "error.hpp"


SparseMtx::SparseMtx(int n, int m, double* et, int* cn, int* da)
{
  nrows = n;
  length = m;
  sra = new double [length];
  clm = new int [length];
  fnz = new int [nrows +1];
  for (int i = 0; i < length; i++)
  {
    sra[i] = et[i];
    clm[i] = cn[i];
  }
  for (int i = 0; i <= nrows; i++)  fnz[i] = da[i];
}

SparseMtx::SparseMtx(int n, int m)
{
  nrows = n;
  length = m;
  sra = new double [length];
  clm = new int [length];
  fnz = new int [nrows +1];

  for (int i =0; i< length; i++)
  {
    sra[i] = 0;
    clm[i] = 0;
  }
  for (int i =0; i <= nrows; i++)  fnz[i] = 0;
}

SparseMtx::SparseMtx(const SparseMtx& mat)  // Copy Constructor
{
  nrows = mat.nrows;
  length = mat.length;
  sra = new double [length];
  clm = new int [length];
  fnz = new int [nrows +1];
  for (int i = 0; i < length; i++) {
    sra[i] = mat[i];
    clm[i] = mat.clm[i];
  }
  for (int i =  0; i <= nrows; i++) fnz[i] = mat.fnz[i];
}


SparseMtx& SparseMtx::operator=(const SparseMtx & ssm)
{
  if(nrows != ssm.nrows || length != ssm.length)
  error("Bad matrix sizes in SparseMtx::operator=()");
  for (int i = 0; i < length; i++)
  {
    sra[i]  = ssm[i];
    clm[i]  = ssm.clm[i];
  }
  for (int i = 0; i <= nrows; i++) fnz[i] = ssm.fnz[i];
  return *this;
}

Vector SparseMtx::operator*(const Vector& vec) const
{
  if (nrows != vec.size())
  error("Matrix-Vector sizes do not match in SparseMtx::operator*(). ");
  Vector tm(nrows);
  for (int i = 0; i < nrows; i++)
  for (int j = fnz[i]; j < fnz[i +1]; j++)
  tm[i] += sra[j]*vec[clm[j]];
  return tm;
}
