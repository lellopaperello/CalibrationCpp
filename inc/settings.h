// Header for the Settings Class.

#ifndef SETTINGS_H
#define SETTINGS_H


using namespace std;

// Structures
struct phi_t{
  string  name;
  int     N;
  int     K;
  double  vec[10]; // Just for hardcoded input
};

class Settings {
public:
  // Project name
  string name = "InsertNameHere";

  // Data Set
  string type = "InsertTypeHere";
  testCase_t testCase;

  // Data Analysis
  string         model = "InsertModelHere";
  phi_t          phi;
  vector<double> pi;

  // Ambient parameters
  double rho_a = 1.225;
  double mu = 1.715e-5;
  double g = 9.81;
private:
protected:
};

#endif
