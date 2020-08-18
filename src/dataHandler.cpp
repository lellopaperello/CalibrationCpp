// Class implementation.
// Contains the routines to generate the data set (starting from a Test Case or
// a relation present in the Literature) and to organize it.

#include "dataHandler.h"

data_t DataHandler::GenerateTestCase(const testCase_t testCase,
                                     const vector<double>::size_type& nData,
                                     const vector<double>& D) {

testName_t testName = GetTestCase(testCase.name);
Literature lit(testCase.model);

vector<vector<double>> phi (nData, vector<double> (testCase.nParam));
default_random_engine generator;

double sigmaD, sigmaVt;

data_t data;
data.N = nData * D.size();
data.dv.resize(data.N);
data.vt.resize(data.N);
data.sigma.resize(data.N);

  switch (testName) {

    case MONOMODAL: {
      normal_distribution<double> monomodal(testCase.mu, testCase.sigma); // Only 1 mode allowed

      // Generate the random PHI samples
      for (vector<double>::size_type i = 0; i < nData; i++) {
        for (vector<double>::size_type j = 0; j < testCase.nParam; j++) {
          phi[i][j] = monomodal(generator);
        }
      }

      // Measurement error on the Terminal Velocity
      if (testCase.sigma == 0) {
        sigmaVt = 1.0e-4;
      } else {
        sigmaVt = testCase.sigma / 10;
      }

      // Measurement error on the Diameter
      sigmaD = (D.back() - D.front()) / (0.5 * D.size());

      data.name = "Standard Distribution";
      data.longname = data.name + ", mu = " + to_string(testCase.mu)
                                + ", sigma = " + to_string(testCase.sigma);
    break;
    }

    // No valid model selected
    default:
    // ERROR MESSAGE
    break;
  }

  for (vector<double>::size_type i = 0; i < D.size(); i++) {
    normal_distribution<double> diameter(D[i], sigmaD);

    for (size_t j = 0; j < nData; j++) {
      size_t ind = i * nData + j;

      // Generate uniformly distributed diameters around the given one
      data.dv[ind] = diameter(generator);

      // Generate terminal velocities from the selected model
      data.vt[ind] = lit.CalculateVt(data.dv[ind], phi[j]);

      // Measurement error on the Terminal Velocity
      data.sigma[ind] = sigmaVt;
    }
  }

  return data;
}


// Test Case selector
DataHandler::testName_t DataHandler::GetTestCase(const string& TESTCASE) {
  if ((TESTCASE == "MONOMODAL") | (TESTCASE == "Monomodal") |
      (TESTCASE == "monomodal") | (TESTCASE == "1") | (TESTCASE == "Test 1")) {
    return MONOMODAL;
  }
  else if (TESTCASE == "BIMODAL") {
    return BIMODAL;
  }
  else {
    cout << "Invalid model selected."; // CAPIRE COME USARE GLI ERRORI
    return MONOMODAL;
  }

}
