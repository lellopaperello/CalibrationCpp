// Class implementation.
// Contains the routines to generate the data set (starting from a Test Case or
// a relation present in the Literature) and to organize it.

#include "dataHandler.h"

data_t DataHandler::GenerateTestCase(const testCase_t testCase,
                                     const vector<double>::size_type& nData,
                                     const vector<double>& D) {

testName_t testName = GetTestCase(testCase.name);

vector<vector<double>>   phi (nData, vector<double> (testCase.nParam));
function<double(double)> vtMean;
double                   vtSigma;
default_random_engine    generator;

double sigmaD, sigmaVt;

data_t data;
data.N = nData * D.size();
data.dv.resize(data.N);
data.vt.resize(data.N);
data.sigma.resize(data.N);

  switch (testName) {

    case MONOMODAL: {
      vector<normal_distribution<double>> monomodal;
      for (vector<vector<double>>::size_type j = 0; j < testCase.nParam; j++) {
        monomodal.push_back(normal_distribution<double>
                           (testCase.mu[j][0], testCase.sigma[j][0]));
      }

      // Generate the random PHI samples
      for (vector<double>::size_type i = 0; i < nData; i++) {
        for (vector<double>::size_type j = 0; j < testCase.nParam; j++) {
          phi[i][j] = monomodal[j](generator);
        }
      }

      // Measurement error on the Terminal Velocity
      if (testCase.sigma[0][0] == 0) {
        sigmaVt = 1.0e-4;
      } else {
        sigmaVt = testCase.sigma[0][0];
      }

      // Measurement error on the Diameter
      sigmaD = (D.back() - D.front()) / (0.5 * D.size());

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

      // Measurement error on the Terminal Velocity
      sigmaVt = 0.01;

      // Measurement error on the Diameter
      sigmaD = (D.back() - D.front()) / (0.5 * D.size());

      data.name = "Brandes -1°C";
      data.longname = "Brandes - 2008, T = -1°C";

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
      if (testName > 10) {
        Literature lit(testCase.model);
        data.vt[ind] = lit.CalculateVt(data.dv[ind], phi[j]);
      } else {
        normal_distribution<double> velocity (vtMean(data.dv[ind]), vtSigma);
        data.vt[ind] = velocity(generator);
      }

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
