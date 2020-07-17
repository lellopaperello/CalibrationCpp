// Header for the Literature Class.

#ifndef LITERATURE_H
#define LITERATURE_H

  #include <math.h>
  #include <functional>

  #include "io.h"
  #include "dataHandler.h"
  
  using namespace std;

class Literature
{
public:
  // Constructors
  Literature(const string& MODEL);
  Literature(const string& MODEL, double& RHO_A, double& MU, double& G);

  // Destructor
  ~Literature() {}

  // Methods
  double CalculateVt(double dv, const vector<double>& phi);

protected:
private:
  // Enum
  enum model_t { GANSER, HOLTZERSOMMERFELD };

  // Attributes
  model_t model;  // chosen model     [ string ]
  double rho_a;   // air density      [ kg/m^3 ]
  double mu;      // air viscosity    [ kg/ms ]
  double g;       // gravity field    [ N/kg ]

  // Methods
  double rho_snow(double& D);
  model_t GetModel(const string& MODEL);

};

#endif
