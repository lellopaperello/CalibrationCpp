// Header for the Overloading of the required operators of the vector class.

#ifndef VECTOR_OVERLOAD_H
#define VECTOR_OVERLOAD_H
  #include <vector>

  using namespace std;

// Vector<double>
vector<long double> operator+(const vector<long double>& v1,
                              const vector<long double>& v2);

void operator+=(vector<long double>& lhs, const vector<long double>& rhs);

vector<long double> operator*(const vector<long double>& v1, double& alpha);

#endif // VECTOR_OVERLOAD_H
