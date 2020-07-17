// Class implementation.
// Contains the routines to generate the data set (starting from a Test Case or
// a relation present in the Literature) and to organize it.

data_t DataHandler::GenerateTestCase() {

  switch (testCase) {
    case /* value */: {
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
