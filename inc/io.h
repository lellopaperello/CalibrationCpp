// Header for the Input/Output Class.

#ifndef IO_H
#define IO_H

  #include <iostream>
  #include <vector>
  #include <string>

  using namespace std;

// Enum
enum model_t { GANSER };

// Structures
struct data_t {
  string name;
  string longname;
  double N;
  double* vt;
  double* dv;
  double* sigma;
};

struct phi_t{
  string name;
  double* vec;
};

class io
{
public:
  // Meethods
  void printvec(const vector<int>& v);
};
#endif // IO_H
