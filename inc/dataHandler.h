testCase// Header for the DataHandler Class.

#ifndef DATA_HANDLER_H
#define DATA_HANDLER_H

#include <random>
#include <string>

#include "io.h"
#include "literature.h"

using namespace std;

// Structures
struct data_t {
  string  name;
  string  longname;
  double  N;
  double* vt;
  double* dv;
  double* sigma;
};


class DataHandler {
public:
  data_t GenerateTestCase();

private:
protected:
  // Enum
  enum testCase_t { /* Put Names here */ };

  // Methods
  model_t GetTestCase(const string& TESTCASE);

};

#endif
