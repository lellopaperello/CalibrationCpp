// Extension of the program that manages the input and output.

#include "io.h"

void IO::printvec(const vector<int>& v) {
  for (auto i : v) {
    cout << i;
  }
  cout << '\n';
}

void IO::printPosterior(const vector<double>& v, const string& file) {
  ofstream out;
  out.open(file);

  for (auto i : v) {
    out << i << '\n';
  }

  out.close();
}
