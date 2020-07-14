// Header for the Literature Class.

#ifndef LITERATURE_H
#define LITERATURE_H

  #include <math.h>
  #include <functional>
  
  #include "io.h"

  using namespace std;

// Constants
#define rho_a 1.225         // air density               [ kg/m^3 ]
#define mu    1.715e-5      // air viscosity             [ kg/ms ]
#define g     9.81          // gravity field             [ N/kg ]


class Literature
{
public:
  // Methods
  double CalculateVt(double dv, const vector<double>& phi, const model_t& model);

protected:
private:
};

#endif
