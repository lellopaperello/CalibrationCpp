// Extension of the program that manages the input and output.

#include "io.h"

Settings loadSettings(string configFile) {
  Config cfg;
  Settings settings;

  // Read the file. If there is an error, report it and exit.
  try
  {
    cfg.readFile(configFile);
  }
  catch(const FileIOException &fioex)
  {
    cerr << "I/O error while reading file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const ParseException &pex)
  {
    cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
         << " - " << pex.getError() << endl;
    return(EXIT_FAILURE);
  }

  // Get the Project name.
  try
  {
    settings.name = cfg.lookup("name");
    cout << "Store name: " << settings.name << endl << endl;
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No 'name' setting in configuration file." << endl;
  }

  const Setting &root = cfg.getRoot();

  // Get the Data Set information.
  try
  {
    const Setting &dataSet = root["dataSet"];

    // Data Set type
    try
    {
      settings.type = dataSet.lookup("type")
    }












}


void IO::printvec(const vector<int>& v) {
  for (auto i : v) {
    cout << i;
  }
  cout << '\n';
}


void IO::printPosterior(const vector<double>& v, const string& file) {
  ofstream out;
  out.open(file);

  for (auto i : v) {
    out << i << '\n';
  }

  out.close();
}
