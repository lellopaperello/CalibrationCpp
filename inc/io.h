// Header for the Input/Output Class.

#ifndef IO_H
#define IO_H

  #include <iostream>
  #include <fstream>
  #include <vector>
  #include <string>
  #include <iomanip>
  #include <cstdlib>
  #include <libconfig.h++>

  #include "settings.h"

  using namespace std;
  using namespace libconfig;

class IO
{
public:
  // Methods
  Settings loadSettings(const char *configFile);
  void printProgress(int i, int &prog);
  void printPosterior(const vector<double>& v, const string& file);

  static void printVec(const vector<int>& v);
  static void printVec(const vector<double>& v);
  static void printVec(const vector<float>& v);
};

#endif // IO_H
