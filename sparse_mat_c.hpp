#ifndef SPARSE_MAT_C_H
#define SPARSE_MAT_C_H
using namespace std;
#include "Vector.hpp"
#include "abstract_mat_c.hpp"

// Compressed Sparse Row Format
class SparseMtx: public AbsMtx
{
  private:
  int length;             // Number of non-zero entries of the original matrix
  double* sra;           // Array for storing the non-zero entries
  int* clm;              // Column indexes in matrix of the entries in sra
  int* fnz;              // Position in sra of first non-zero entries of each row

  public:
  SparseMtx(int n, int m, double* t, int* c, int* f);
  // n: number of rows (and columns) of the original matrix
  // m: length of array sra for nonzero entries.
  // t: nonzero entries of the original matrix
  // c: colunm indexes (in the original matrix) of entries in sra
  // f: index in sra of first nonzero entry in each row

  SparseMtx(int n, int m);                  // Initialize all entris to Zero
  // n: number of rows and columns
  // m: number of nonzero entries

  SparseMtx(const SparseMtx&);                   // Copy constructor
  ~SparseMtx(){ delete[] sra; delete[] fnz; delete[] clm; } // Destructor
  SparseMtx& operator=(const SparseMtx&);        // Overload of "="-Operator
  Vector operator*(const Vector&) const;         // Matrix-Vector Product
  double& operator[](int i) const { return sra[i]; }  // Subscripting
  int& getfnz(int i) const { return fnz[i]; } // First !0 entry of each row
  int& getclm(int i) const { return clm[i]; }   // Column Index
};

#endif
