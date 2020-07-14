// Class implementation.
// Contains various models and empirical realtions from the Literature.

#include "literature.h"
#include "solvers.h"


double Literature::CalculateVt(double dv, const vector<double>& phi,
                               const model_t& model) {

  // Number of shape parameters

  auto N = phi.size();
  double Ar;
  function<double(double)> equilibrium;
  function<double(double)> cD;

  // Implementation of different models
  if (model == GANSER) {
        double Phi, dn, K1, K2, Re_v, RE_v;

        // Indipendent shape parameters for Ganser's model
        Phi = phi[0];     dn = phi[1] * dv;

        // Aspect ratio ~ V / S
        Ar = 2/3 * pow(dv, 3) / pow(dn, 2);

        // Stokes shape factor
        K1 = pow(dn/(3*dv) + 2/(3*sqrt(Phi)), -1);

        // Newton shape factor
        K2 = pow(10, (1.8148*(pow(log(Phi), 0.5743))));

        // Reynolds number per unit velocity
        Re_v = rho_a * dv / mu;

        // Generalized Reynolds number (per unit velocity)
        RE_v = Re_v * K1 * K2;

        // Drag coefficient model
        cD = [&](double v)
        {
          return K2* (24 / (RE_v*v) * (1 + 0.1118 * pow(RE_v*v, 0.6567))
                   + 0.4305 / (1 + 3305 / (RE_v*v)));
        };
  }

  // Equilibrium equation
  equilibrium = [&] (double vt)
  {
    double rho_snow;
    return  1/2 * pow(vt, 2) * cD(vt) + (rho_a - rho_snow)/rho_a * Ar * g;

  };

  // Solution
  BisectionSolver _solver(equilibrium, 1.0e-10, 1.0e3);
  return _solver.solve();

}
