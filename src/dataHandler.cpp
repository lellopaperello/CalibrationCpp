// Class implementation.
// Contains the routines to generate the data set (starting from a Test Case or
// a relation present in the Literature) and to organize it.

data_t DataHandler::GenerateTestCase() {

vactor<vector<double>> phi (nParam, vector<double> (nData));

  switch (testCase) {
    // Generate the random PHI samples
    case MONOMODAL: {
      default_random_engine generator;
      normal_distribution<double> monomodal(mu, sigma);

      for (vector<double>::size_type i = 0; i < nParam; i++) {
        for (vector<double>::size_type j = 0; j < nData; j++) {
          phi[i][j] = monomodal(generator);
        }
      }

      name = "Standard Distribution"
      longname = name + ", mu = " + to_string(mu)
                      + ", sigma = " + to_string(sigma);
    break;
    }

    // No valid model selected
    default:
    // ERROR MESSAGE
    break;
  }









}


// Test Case selector
DataHandler::testCase_t DataHandler::GetModel(const string& TESTCASE) {
  if ((TESTCASE == /* Insert name here */) | (TESTCASE == /* Insert name here */)) {
    return /* Insert name here */;
  }
  else if (TESTCASE == /* Insert name here */) {
    return /* Insert name here */;
  }
  else {
    throw "Invalid model selected."; // CAPIRE COME USARE GLI ERRORI
  }

}
