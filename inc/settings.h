// Header for the Settings Class.

#ifndef SETTINGS_H
#define SETTINGS_H


using namespace std;

// Structures
struct phi_t{
  string          name;
  string          latex;
  int             K;
  vector<double>  vec;
};

struct data_t {
  string         name;
  string         longname;
  int            N;
  vector<double> vt;
  vector<double> dv;
  vector<double> sigma;
};

struct testCase_t {
  string                 name;
  string                 model;
  int                    nParam;
  int                    nData;
  vector<double>         D;
  vector<double>         pi;
  vector<vector<double>> mu;      // 1st index: Parameter
  vector<vector<double>> sigma;   // 2nd index: Mode
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
  vector<phi_t>  phi;
  vector<double> pi;

  // Ambient parameters
  double rho_a = 1.225;
  double mu = 1.715e-5;
  double g = 9.81;
private:
protected:
};

#endif
