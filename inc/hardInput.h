// Hard-coded input for testing purpose.

#ifndef HARD_INPUT_H
#define HARD_INPUT_H

  #include "dataHandler.h"

  using namespace std;

testCase_t testCase = {
  .name = "MONOMODAL",
  .model = "GANSER",
  .nParam = 2,
  .mu = 0.5,
  .sigma = 0.01
};

vector<double>::size_type nData = 10;

vector<double> D = {1e-3, 1.25e-3, 1.5e-3};


vector<phi_t> phi = {
  {
    .name = "Phi",
    .N = 10,
    .vec = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
  }, {
    .name = "dn / dv",
    .N = 10,
    .vec = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
  }
};

string model = "GANSER";

#endif // HARD_INPUT_H
