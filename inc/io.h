// Header for the Input/Output Class.

#ifndef IO_H
#define IO_H

  #include <iostream>
  #include <fstream>
  #include <vector>
  #include <string>

  using namespace std;

class IO
{
public:
  // Methods
  Settings loadSettings(string configFile);
  void printvec(const vector<int>& v);
  void printPosterior(const vector<double>& v, const string& file);
};

#endif // IO_H
