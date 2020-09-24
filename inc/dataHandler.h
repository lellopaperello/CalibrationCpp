// Header for the DataHandler Class.

#ifndef DATA_HANDLER_H
#define DATA_HANDLER_H

#include <random>
#include <string>
#include <vector>

#include "io.h"
#include "literature.h"

using namespace std;

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
  // Enum: Use identifiers > 10 for test cases taken from the literature
  enum testName_t { MONOMODAL,
                    BIMODAL,
                    BRANDES1 = 11 };

  // Methods
  testName_t GetTestCase(const string& TESTCASE);

};

#endif
