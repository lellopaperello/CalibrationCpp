// Calibration of snowflakes shape parameters using Bayesian approach.

#include "posterior.h"
#include "dataHandler.h"
#include "io.h"

int main(int argc, char const *argv[]) {

  DataHandler dataHandler;
  IO          io;

  // Loading Input Settings
  Settings settings = io.loadSettings("Brandes.cfg");

  string         postOutFIle = "res/brandes.txt";
  // Generation of an artificial dataset ---------------------------------------
  data_t data = dataHandler.GenerateTestCase(settings.testCase);


  // Subdivision of the data set -----------------------------------------------

  /*   TO BE IMPLEMENTED   */


  // Data Analysis - Posterior generation --------------------------------------
  Posterior posterior(data, settings.phi, settings.pi, settings.model);
  vector<double> post = posterior.bruteForce();


  // Saving results
  io.printPosterior(post, postOutFIle);

  return 0;
}
