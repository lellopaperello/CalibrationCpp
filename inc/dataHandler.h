// Header for the DataHandler Class.

#ifndef DATA_HANDLER_H
#define DATA_HANDLER_H

#include <random>
#include <string>
#include <vector>

#include "io.h"
#include "literature.h"

using namespace std;

// Structures
struct data_t {
  string                     name;
  string                     longname;
  vector<double> ::size_type N;
  vector<double>             vt;
  vector<double>             dv;
  vector<double>             sigma;
};

struct testCase_t {
  string                     name;
  string                     model;
  vector<double> ::size_type nParam;
  double                     mu;
  double                     sigma;
};

class DataHandler {
public:
  // Constructor
  DataHandler() {}

  // Destructor
  ~DataHandler() {}

  // Methods
  data_t GenerateTestCase(const testCase_t testCase,
                          const vector<double> ::size_type& nData,
                          const vector<double>& D);

private:
protected:
  // Enum
  enum testName_t { MONOMODAL, BIMODAL };

  // Methods
  testName_t GetTestCase(const string& TESTCASE);

};

#endif
