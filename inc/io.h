// Header for the Input/Output Class.

#ifndef IO_H
#define IO_H

  #include <iostream>
  #include <fstream>
  #include <vector>
  #include <string>

  using namespace std;

// Structures
struct phi_t{
  string  name;
  int     N;
  int     K;
  double  vec[10]; // Just for hardcoded input
};

class IO
{
public:
  // Methods
  void printvec(const vector<int>& v);
  void printPosterior(const vector<double>& v, const string& file);
};

#endif // IO_H
