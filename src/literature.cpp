// Class implementation.
// Contains various models and empirical realtions from the Literature.

#include "literature.h"

// Default constructor: Set only the model
Literature::Literature(const string& MODEL) {
  model = GetModel(MODEL);
  rho_a = 1.225;
  mu    = 1.715e-5;
  g     = 9.81;
};

// Constructor with parameters: Set the model and the ambient parameters
Literature::Literature(const string& MODEL, double& RHO_A, double& MU, double& G) {
  model = GetModel(MODEL);
  rho_a = RHO_A;
  mu    = MU;
  g     = G;
};

// Model selector
Literature::model_t Literature::GetModel(const string& MODEL) {
  if ((MODEL == "Ganser") | (MODEL == "ganser") | (MODEL == "GANSER")) {
    return GANSER;
  }
  else if (MODEL == "H&S") {
    return HOLTZERSOMMERFELD;
  }
  else {
    cout << "Invalid model selected."; // CAPIRE COME USARE GLI ERRORI
    return GANSER;
  }
}

// Terminal Velocity calculator
double Literature::CalculateVt(const double dv, const vector<double>& phi) {

  // Number of shape parameters

  auto N = phi.size();
  double Ar;
  function<double(double)> equilibrium;
  function<double(double)> cD;

  // Implementation of different models
  switch (model) {

    case GANSER: {
      if (phi.size() != 2) {
        cout << "Ganser model requirest 2 parameters!" << '\n';
      }
      double Phi, dn, K1, K2, Re_v, RE_v;

      // Indipendent shape parameters for Ganser's model
      Phi = phi[0];     dn = phi[1] * dv;

      // Aspect ratio ~ V / S
      Ar = 2/(double)3 * pow(dv, 3) / pow(dn, 2);

      // Stokes shape factor
      K1 = pow(dn/(3*dv) + 2/(3*sqrt(Phi)), -1);

      // Newton shape factor
      K2 = pow(10, (1.8148*(pow(abs(log(Phi)), 0.5743))));

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

      break;
    }

    // No valid model selected
    default:
    // ERROR MESSAGE
    break;
  }

  // Equilibrium equation
  equilibrium = [&] (double vt)
  {
    return  1/(double)2 * pow(vt, 2) * cD(vt) + (rho_a - rho_snow(dv))/rho_a * Ar * g;

  };

  // Solution
  BisectionSolver _solver(equilibrium, 1.0e-10, 1.0e3);
  return _solver.solve();

}


double Literature::rho_snow(const double& D) {
  // Brandes model for the calculation of the snowflakes density as function
  // of the ground temperature ~ Brandes et al. 2007,
  // "A Statistical and Physical Description of Hydrometeor Distributions
  //  in Colorado Snowstorms Using a Video Disdrometer"

  // D: equivalent volume diameter [ m ]
  // rho: snowflake density [kg / m^3]

  double D_mm = D * 1.0e3;
  return 0.178e3 * pow(D_mm, -0.922);
}
