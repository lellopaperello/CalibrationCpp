// Calibration of snowflakes shape parameters using Bayesian approach.

#include "posterior.h"
#include "dataHandler.h"
#include "io.h"

int main(int argc, char const *argv[]) {

  DataHandler dataHandler;
  IO          io;

  // Loading Input Settings
  Settings settings = io.loadSettings("config/gaTest.cfg");

  // string postOutFIle = "res/validation.txt";
  // Generation / Loading of an artificial dataset -----------------------------
  data_t data = dataHandler.GenerateTestCase(settings.testCase);
  // data_t data = dataHandler.LoadData("config/dataset.dat", false);


  // Data Analysis - Posterior generation --------------------------------------
  Posterior posterior(data, settings.phi, settings.pi, settings.model);
  // vector<double> post = posterior.bruteForce();
  posterior.GeneticAlgorithm();

  // Saving results
  // io.printPosterior(post, postOutFIle);

  return 0;
}
