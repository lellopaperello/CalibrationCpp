// Class implementation.
// Contains various implementations of solvers for the approximate solution of
// non-linear equations.

#include "solvers.h"

double BisectionSolver::solve() {
	double medium = (_right + _left) / 2;
	if (_func(medium) == 0)
		return medium;
	while (true) {
		medium = (_right + _left) / 2;
		if (medium == _right || medium == _left)
			break;
		if (_func(medium) * _func(_right) < 0)
			_left = medium;
		else
			_right = medium;
	}
	return medium;
}
std::string BisectionSolver::name() {
	return "Bisection method";
}
