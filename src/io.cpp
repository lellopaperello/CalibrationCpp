// Extension of the program that manages the input and output.

#include "io.h"

void io::printvec(const vector<int>& v) {
  for (auto i : v) {
    cout << i;
  }
  cout << '\n';
}
