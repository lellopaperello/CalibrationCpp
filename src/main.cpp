// Calibration of snowflakes shape parameters using Bayesian approach.

#include "posterior.h"
#include "dataHandler.h"
#include "io.h"

int main(int argc, char const *argv[]) {

  DataHandler dataHandler;
  Settings settings;
  IO          io;

  // Loading Input Settings
  if (argc == 2) {
    settings = io.loadSettings(argv[1]);
  } else if (argc == 1) {
    settings = io.loadSettings("config/config.cfg");
  } else { // argc > 2
    cout << "Error: too many input arguments.";
    return 0;
  }


  // Generation / Loading of an artificial dataset -----------------------------
  data_t data;
  if (settings.type == "generate") {
    data = dataHandler.GenerateTestCase(settings.testCase);
  } else if (settings.type == "load") {
    data = dataHandler.LoadData(settings.dataInFile, settings.printData);
  }


  // Data Analysis - Posterior generation --------------------------------------
  Posterior posterior(data, settings.phi, settings.K, settings.pi, settings.model);
  if (settings.approach == "bruteforce") {
    vector<double> post = posterior.bruteForce();
    io.printPosterior(post, settings.postOutFile);
  } else if (settings.approach == "GA") {
    posterior.GeneticAlgorithm(settings.gaInputFile, settings.gaOutputFile);
  }

  return 0;
}
