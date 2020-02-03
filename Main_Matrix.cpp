#include <iostream>
using namespace std;
#include "Vector.hpp"
#include "full_mat_c.hpp"
#include "sparse_mat_c.hpp"
#include "abstract_mat_c.hpp"

int main()
{
  cout << " ******** Testing for Full matrix  ******** \n";
  int n = 2, m = 2;
  FullMtx A(n,m);
  /*for (int i = 0; i < n; i++)
  for (int j = 0; j < m; j++)
  A[i][j] =  n/(abs(i - j) + 1.0);*/
  A = {{4,1},{1,3}};
  cout << "**********" << endl;
  cout <<"Matrix A is: "<< endl << A << endl;
  cout << "**********" << endl;
  FullMtx B(A);
  B[0][1] = -1;
  cout << "Matrix B is " << endl << B << endl;
  cout << "**********" << endl;
  FullMtx C(A*B);
  cout << "Matrix C is " << endl << C << endl;
  cout << "**********" << endl;
  FullMtx D(A+B);
  cout << "Matrix D is " << endl << D << endl;
  cout << "**********" << endl;
  Vector v1(n,0.0);
  double alpha = 90.90909090;
  cout << "Alpha is : " << alpha << endl;
  v1[0] = alpha/1000;
  v1[1] = 7*alpha/1000;
  cout <<"Vector 1: "<< endl <<  v1 << endl;
  Vector VTS = v1;
  cout << "Vector VTS: "<< endl << VTS << endl;

  cout << " ******** Testing for Abs matrix  ******** \n";
  int prec = 0;
  unsigned int iter = 2000;
  double eps = 1.0e-14;
  Vector vec2(n,0.0);
  vec2[0] = (2);
  vec2[1] = (1);
  cout << "Vector 2: "<< endl << vec2 << endl;
  Vector b(n,1);
  b[0] = 1; b[1]=2;
  //cout << (A*v1).size()<<endl;
  //cout << v1.size()<<endl;
  //cout << A.getnrows()<<endl;
  int ret =  A.CG(vec2, b, eps,iter,prec); // Ax = b <=> x = A^(-1)b
  if (ret == 0) cout << "CG returned successfully " << endl;
  cout << "CGsolution =  "<< endl << vec2 << endl;
  string iteration = " ";
  (iter == 1) ? iteration = " iteration " : iteration = " iterations ";
  cout << iter << iteration << " used; " << endl;
  cout << "Residual in CG = " << eps << " ; "<< endl;
  cout << "True Solution is : " << endl << VTS << endl;
  cout << "CG Solution is : " << endl << vec2 << endl;
  cout << "Difference between true solution & CG solution is : " <<
  endl << vec2-VTS<< endl;
  cout << "Ture error in CG = " << (vec2-v1).maxnorm() << " . " << endl;
  cout << "Vector b is: "<< endl<<  b << endl;
  cout << "Matrix-Vector Product between A and x_conjuguate Gradient is: "<<
  endl << A*vec2<< endl;


  cout << "******** Testing for Sparse Matrix  ******** \n";
  SparseMtx sm1(n, (n*n));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      sm1[i*n + j] = A[i][j];
      sm1.getclm(i*n + j) = j;
    }
    sm1.getfnz(i) = i*n;
  }
  sm1.getfnz(n) = n*n;
  iter = n;
  eps = 1.0e-14;
  //vec2.reset();
  ret =  sm1.CG(vec2, sm1*v1,eps,iter,prec);
  if (ret == 0) cout << "CG returned successfully\n";
  // cout << "True solution is: " << v1 << " ";
  // cout << "CG solution =  " << vec2 << "\n";
  cout << iter << " iterations used. " ;
  cout << "Residual in CG = " << eps << "  " ;
  cout << "Ture error in CG = " << (vec2-v1).maxnorm() << "\n";

}
