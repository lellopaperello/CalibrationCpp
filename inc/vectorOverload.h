// Header for the Overloading of the required operators of the vector class.

#ifndef VECTOR_OVERLOAD_H
#define VECTOR_OVERLOAD_H
  #include <vector>

  using namespace std;

// Vector<double>
vector<double> operator+(const vector<double>& v1, const vector<double>& v2) {
  // Overloading of the operator+ for the vector<double> class.

  // Vectors must be of the same size!          ERROR HANDLING!

  vector<double>::size_type N = v1.size();
  vector<double> result (N);

  for (vector<double>::size_type i = 0; i < N; i++) {
    result[i] = v1[i] + v2[i];
  }

  return result;
}

vector<double> operator-(const vector<double>& v1, const vector<double>& v2) {
  // Overloading of the operator- for the vector<double> class.

  // Vectors must be of the same size!          ERROR HANDLING!

  vector<double>::size_type N = v1.size();
  vector<double> result (N);

  for (vector<double>::size_type i = 0; i < N; i++) {
    result[i] = v1[i] - v2[i];
  }

  return result;
}

vector<double> operator*(const vector<double>& v1, double& alpha) {
  // Overloading of the operator* for the vector<double> class.
  // Product vector by scalar.

  // Vectors must be of the same size!          ERROR HANDLING!

  vector<double>::size_type N = v1.size();
  vector<double> result (N);

  for (vector<double>::size_type i = 0; i < N; i++) {
    result[i] = v1[i] * alpha;
  }

  return result;
}

vector<double> operator*(const vector<double>& v1, const vector<double>& v2) {
  // Overloading of the operator* for the vector<double> class.
  // Product vector by vector, element by element.

  // Vectors must be of the same size!          ERROR HANDLING!

  vector<double>::size_type N = v1.size();
  vector<double> result (N);

  for (vector<double>::size_type i = 0; i < N; i++) {
    result[i] = v1[i] * v2[i];
  }

  return result;
}

vector<double> operator/(const vector<double>& v1, double& alpha) {
  // Overloading of the operator/ for the vector<double> class.
  // Division vector by scalar.

  // Vectors must be of the same size!          ERROR HANDLING!

  vector<double>::size_type N = v1.size();
  vector<double> result (N);

  for (vector<double>::size_type i = 0; i < N; i++) {
    result[i] = v1[i] * alpha;
  }

  return result;
}

vector<double> operator/(const vector<double>& v1, const vector<double>& v2) {
  // Overloading of the operator/ for the vector<double> class.
  // Division vector by vector, element by element.

  // Vectors must be of the same size!          ERROR HANDLING!

  vector<double>::size_type N = v1.size();
  vector<double> result (N);

  for (vector<double>::size_type i = 0; i < N; i++) {
    result[i] = v1[i] * v2[i];
  }

  return result;
}


#endif // VECTOR_OVERLOAD_H
