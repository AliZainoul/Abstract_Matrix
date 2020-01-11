#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
#include <cassert>
using namespace std;

class Vector
{
 private:
  double* M_data;
  int M_size;

 public:
  Vector(); // Default Constructor
  // Vector(int i); // Constructor 1
  Vector(int i, double v=0); // Constructor 2
  Vector(initializer_list<double> lst); // Constructor 3
  Vector(Vector const &v); // Copy Constructor
  Vector& operator=(Vector const &v); // Assignment Operator
  // Subscript Operator lvalue (for left value)
  double& operator[](int i){return M_data[i];}
  // Subscript Operator const
  double operator[](int i) const{return M_data[i];}
  Vector(Vector&& v); // Move constructor
  Vector& operator=(Vector&& v); // Move Assignment Operator
  friend std::ostream& operator << (std::ostream&, const Vector&); //'<<'Overload
  void reset(double = 0);
 ~Vector(); // Destructor
  void affiche();
  int size() const { return M_size; } // Return Size of the Vector
  double maxnorm() const; // L-infinity norm
  double twonorm() const; // L-2 norm
  Vector& operator+=(const Vector& vec); //
  Vector& operator-=(const Vector& vec);
  friend Vector operator+(const Vector&); // Unary operator, V += W
  friend Vector operator-(const Vector&); // Unary operator, V -= W
  friend Vector operator+(const Vector&, const Vector&);// Binary Operator,V=V1+V2
  friend Vector operator-(const Vector&, const Vector&);// Binary Operator,V=V1-V2
  friend Vector operator*(double, const Vector&); // Scalar-Vector Product
  friend Vector operator*(const Vector&, double); // Vector-Scalar Product
  friend double dot(const Vector&, const Vector&); // Dot product
};



#endif
