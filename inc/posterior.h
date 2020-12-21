// Header for the Posterior Class.

#ifndef POSTERIOR_H
#define POSTERIOR_H

  #include <string>
  #include <vector>
  #include <algorithm>
  #include <limits>
  #include <omp.h>

  #include "io.h"
  #include "dataHandler.h"
  #include "literature.h"
  #include "vectorOverload.h"

  #include "ga/ga.h"
  #include <ga/std_stream.h>
  #include <ga/GARealGenome.h>


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
  void GeneticAlgorithm (const string& gaInputFile = "none",
                         const string& gaOutputFile = "none");

  // Useful functions
  static void ind2sub(const vector<int>& size, int ind,  vector<int>& sub);
  static void sub2ind(const vector<int>& size, const vector<int>& sub,  int& ind);

protected:
  // Structure
  struct userData_t {
    data_t       data;
    int          K;
    vector<int>  Kvec;
    string       model;
    double       base;
  };

private:
  // Attributes
  data_t          data;
  vector<phi_t>   phi;
  vector<double>  pi;
  string          model;

  // Methods
  static vector<int> cumprod(const vector<int>& size);
  static vector<double> normConstraint(vector<double> pi);
  static vector<float> normConstraint(vector<float> pi);
  long double findMax(const vector<vector<long double>>& vec, int s1, int s2);

  static float Objective  (GAGenome &);
  static float Objective2 (GAGenome &);
  static float Objective3 (GAGenome &);
};

#endif // POSTERIOR_H
