// Implementation of the overloading of the required operators of the vector
// class.

#include "vectorOverload.h"

// Vector<double>
vector<long double> operator+(const vector<long double>& v1,
                              const vector<long double>& v2) {
  // Overloading of the operator+ for the vector<long double> class.

  // Vectors must be of the same size!          ERROR HANDLING!

  vector<long double>::size_type N = v1.size();
  vector<long double> result (N);

  for (vector<long double>::size_type i = 0; i < N; i++) {
    result[i] = v1[i] + v2[i];
  }

  return result;
}

void operator+=(vector<long double>& lhs, const vector<long double>& rhs) {
  // Overloading of the operator+= for the vector<long double> class.

  // Vectors must be of the same size!          ERROR HANDLING!

  for (vector<long double>::size_type i = 0; i < rhs.size(); i++) {
    lhs[i] += rhs[i];
  }
}

vector<long double> operator*(const vector<long double>& v1, double& alpha) {
  // Overloading of the operator* for the vector<long double> class.
  // Product vector by scalar.

  // Vectors must be of the same size!          ERROR HANDLING!

  vector<long double>::size_type N = v1.size();
  vector<long double> result (N);

  for (vector<long double>::size_type i = 0; i < N; i++) {
    result[i] = v1[i] * alpha;
  }

  return result;
}
