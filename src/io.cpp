// Extension of the program that manages the input and output.

#include <iostream>
#include <vector>

using namespace std;

void printvec(const vector<int>& v) {
  for (auto i : v) {
    cout << i;
  }
  cout << '\n';
}
