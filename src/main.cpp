// Calibration of snowflakes shape parameters using Bayesian approach.

#include "posterior.h"
#include "dataHandler.h"
#include "io.h"
#include "hardInput.h"


int main(int argc, char const *argv[]) {

  DataHandler dataHandler;
  IO          io;

  vector<double> post;
  string         postOutFIle = "../res/posterior.txt";
  // Generation of an artificial dataset ---------------------------------------
  data_t data = dataHandler.GenerateTestCase(testCase, nData, D);


  // Subdivision of the data set -----------------------------------------------

  /*   TO BE IMPLEMENTED   */


  // Data Analysis - Posterior generation --------------------------------------
  Posterior   posterior(data, phi, model);
  post = posterior.bruteForce();


  // Saving results
  io.printPosterior(post, postOutFIle);

  return 0;
}
