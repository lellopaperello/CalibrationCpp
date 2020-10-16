// Header for the Posterior Class.

#ifndef POSTERIOR_H
#define POSTERIOR_H

  #include <string>
  #include <vector>
  #include <algorithm>
  #include <limits>

  #include "io.h"
  #include "dataHandler.h"
  #include "literature.h"
  #include "vectorOverload.h"

  using namespace std;

class Posterior
{
public:
  // Constructor
  Posterior(const data_t& DATA, const vector<phi_t>& PHI,
            const vector<double>& PI, const string& MODEL);

  // Destructor
  ~Posterior() {}

  // Methods
  vector<double> bruteForce ();

  // Useful functions
  void ind2sub(const vector<int>& size, int ind,  vector<int>& sub);
  void sub2ind(const vector<int>& size, const vector<int>& sub,  int& ind);

protected:
private:
  // Attributes
  data_t          data;
  vector<phi_t>   phi;
  vector<double>  pi;
  string          model;

  // Methods
  vector<int> cumprod(const vector<int>& size);
  vector<double> normConstraint(vector<double> pi);
  long double findMax(const vector<vector<long double>>& vec, int s1, int s2);
};

#endif // POSTERIOR_H
