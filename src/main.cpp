// Calibration of snowflakes shape parameters using Bayesian approach.

#include "posterior.h"
#include "io.h"
#include "solvers.h"

using namespace std;

int main(int argc, char const *argv[]) {

  function<double(double)> f;
  double answer;

  f = [](double x) {
    return x - 5;
  };

  BisectionSolver _solver(f, 0, 10);
  answer = _solver.solve();
  cout << "x = " << answer << '\n';

  return 0;
}
