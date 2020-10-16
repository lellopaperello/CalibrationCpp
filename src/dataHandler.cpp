// Class implementation.
// Contains the routines to generate the data set (starting from a Test Case or
// a relation present in the Literature) and to organize it.

#include "dataHandler.h"

data_t DataHandler::GenerateTestCase(const testCase_t testCase) {

testName_t testName = GetTestCase(testCase.name);

vector<vector<double>>   phi (testCase.nData, vector<double> (testCase.nParam));
function<double(double)> vtMean;
double                   vtSigma;
default_random_engine    generator;

double errD, errVt;

data_t data;
data.N = testCase.nData * testCase.D.size();
data.dv.resize(data.N);
data.vt.resize(data.N);
data.sigma.resize(data.N);

  switch (testName) {

    case MONOMODAL: {
      vector<normal_distribution<double>> monomodal;
      for (int j = 0; j < testCase.nParam; j++) {
        monomodal.push_back(normal_distribution<double>
                           (testCase.mu[j][0], testCase.sigma[j][0]));
      }

      // Generate the random PHI samples
      for (int i = 0; i < testCase.nData; i++) {
        for (int j = 0; j < testCase.nParam; j++) {
          phi[i][j] = monomodal[j](generator);
        }
      }

      // Measurement error on the Terminal Velocity

// FIX!! -----------------------------------------------------------------------

//if (testCase.sigma[0][0] == 0) {
  errVt = 1.0e-4;
//} else {
//  errVt = testCase.sigma[0][0];
//}

// -----------------------------------------------------------------------------

      // Measurement error on the Diameter
      errD = (testCase.D.back() - testCase.D.front())
             / (0.5 * testCase.D.size());

      data.name = "Standard Distribution";
      data.longname = data.name + ", mu = " + to_string(testCase.mu[0][0])
                                + ", sigma = " + to_string(testCase.sigma[0][0]);
      break;
    }

    case BRANDES1:{
      // Function for the mean of the Terminal Velocity
      vtMean = [&] (double dv)
      { // Diameter in [m]
        return 0.87 * pow((dv*1e3), 0.23);
      };

      // Standard deviation of the Terminal Velocity
      vtSigma = 0.22;

      // Measurement error on the Terminal Velocity
      errVt = 0.01;

      // Measurement error on the Diameter
      errD = (testCase.D.back() - testCase.D.front())
             / (0.5 * testCase.D.size());

      data.name = "Brandes -1°C";
      data.longname = "Brandes - 2008, T = -1°C";

      break;
    }

    // No valid model selected
    default:
    // ERROR MESSAGE
    break;
  }

  for (vector<double>::size_type i = 0; i < testCase.D.size(); i++) {
    normal_distribution<double> diameter(testCase.D[i], errD);

    for (int j = 0; j < testCase.nData; j++) {
      int ind = i * testCase.nData + j;

      // Generate uniformly distributed diameters around the given one
      data.dv[ind] = diameter(generator);

      // Generate terminal velocities from the selected model
      if (testName < 10) {
        Literature lit(testCase.model);
        data.vt[ind] = lit.CalculateVt(data.dv[ind], phi[j]);
      } else {
        normal_distribution<double> velocity (vtMean(data.dv[ind]), vtSigma);
        data.vt[ind] = velocity(generator);
      }

      // Measurement error on the Terminal Velocity
      data.sigma[ind] = errVt;
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
  if ((TESTCASE == "BRANDES1") | (TESTCASE == "Brandes1") |
      (TESTCASE == "BRANDES 1") | (TESTCASE == "Brandes 1") |
      (TESTCASE == "brandes1") | (TESTCASE == "brandes 1")) {
    return BRANDES1;
  }
  else if (TESTCASE == "BIMODAL") {
    return BIMODAL;
  }
  else {
    cout << "Invalid model selected."; // CAPIRE COME USARE GLI ERRORI
    return MONOMODAL;
  }

}
