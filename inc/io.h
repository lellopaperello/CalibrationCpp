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
  Settings loadSettings(const char * configFile);
  void printvec(const vector<int>& v);
  void printPosterior(const vector<double>& v, const string& file);
};

#endif // IO_H
