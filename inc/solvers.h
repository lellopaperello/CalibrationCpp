// Header for the Solvers Class and all the implemented solvers.

#ifndef SOLVERS_H
#define SOLVERS_H

#include <string>
#include <functional>

class Solver {
public:
  // Constructor
  explicit Solver(const std::function<double(double)> &func) : _func(func)
  {
  }

  // Destructor
  virtual ~Solver() {}

  // Methods
  virtual double solve() = 0;
  virtual std::string name() = 0;

protected:
  std::function<double(double)> _func;
};

// Child class: Implemented Solvers --------------------------------------------
class BisectionSolver: public Solver{
public:
  // Constructor
  explicit BisectionSolver(const std::function<double(double)> &func, double left, double right)
    : Solver(func), _left(left), _right(right) {
    }

	double solve() override;
	std::string name() override;

private:
  double _left, _right;
};
#endif // SOLVERS_H
