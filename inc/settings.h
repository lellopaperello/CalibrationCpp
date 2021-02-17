// Header for the Settings Class.

#ifndef SETTINGS_H
#define SETTINGS_H


using namespace std;

// Structures
struct phi_t{
  string          name;
  string          latex;
  vector<double>  vec;
  double          step;
};

struct pi_t{
  vector<double>  vec;
  double          step;
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
  testCase_t testCase = {
                          .name = "InsertNameHere",
                          .model = "InsertModelHere",
                          .nParam = 0,
                          .nData = 0,
                          .D = {},
                          .pi = {},
                          .mu = {},
                          .sigma = {}
                        };

  // Data Analysis
  string         approach = "InsertApproachHere";
  string         gaInputFile = "none";
  string         gaOutputFile = "none";
  string         model = "InsertModelHere";
  vector<phi_t>  phi = {};
  int            K = 1;
  pi_t           pi = {};
  double         piStep;

  // Ambient parameters
  double rho_a = 1.225;
  double mu = 1.715e-5;
  double g = 9.81;

  // Input / Output
  string dataInFile  = "dataset.dat";
  bool   printData   = false;
  string outFile     = "screen";
  string postOutFile = "res/posterior.txt";

private:
protected:
};

#endif
