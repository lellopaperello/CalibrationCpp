// Header for the Input/Output Class.

#ifndef IO_H
#define IO_H

  #include <iostream>
  #include <vector>
  #include <string>

  using namespace std;

// Structures
struct phi_t{
  string  name;
  int     N;
  double vec[10]; // Just for hardcoded input
};

class io
{
public:
  // Meethods
  void printvec(const vector<int>& v);
};

#endif // IO_H
