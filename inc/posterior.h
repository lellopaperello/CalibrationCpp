// Header for the Posterior Class.

#ifndef POSTERIOR_H
#define POSTERIOR_H

  #include <string>
  #include <vector>
  #include <algorithm>
  
  #include "io.h"
  #include "dataHandler.h"
  #include "literature.h"

  using namespace std;

class Posterior
{
public:
  // Constructor
  Posterior(const data_t& DATA, const vector<phi_t>& PHI, const string& MODEL);

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
  data_t        data;
  vector<phi_t> phi;
  string        model;

  // Methods
  void cumprod(const vector<int>& size,  vector<int>& k);

};

#endif // POSTERIOR_H
