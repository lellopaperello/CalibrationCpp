// Posterior  = prob(HP|data)
// Likelihood = prob(data|HP)
// prior      = prob(HP)
// -------------------------- Bayes' theorem -----------------------------------
//                   Posterior ~ Likelihood * Prior
// -----------------------------------------------------------------------------
// The Likelihood function is a gaussian difference between the datum and
// the prevision of that datum using the model "model" with the shape
// parameter "phi".

// OUTPUT: - posterior ~ matrix of dimension
//                       ( length(phi(1)) x ... x length(phi(end)) )

// INPUT: - data ~ strucure of array (each element is a different
//                 experimental measurement). must contain:
//               data.dv(i)    = volume equivalent diameter [m]
//               data.vt(i)    = terminal velocity [m/s]
//               data.sigma(i) = std of the datum (experimental uncertainty)

//        - model ~ string (name of the model to be used)

//        - phi ~ array of structure(each element is a different shape
//                parameter for the vt calculation). must contain:
//              phi(i).vec = array of possible values [-]

//        - prior (optional) ~ prior distribution (must be consistent with
//                             the dimensions of the posterior)
//                [default] = uniform
// -----------------------------------------------------------------------------

#include "posterior.h"
#include "literature.h"

// Constructor
Posterior::Posterior(const data_t& DATA, const vector<phi_t>& PHI,
                       const string& MODEL) {
  data = DATA;
  phi = PHI;
  model = MODEL;
};

// Brute Force approach
vector<double> Posterior::bruteForce () {

  // General Declarations

  // Declarations for the Monomodal posterior array [Mono]
  // Storage order: PostMono = [phi1, ..., phiN | nData]
  int            nElMono = 1;               // Number of elements
  vector<int>    sizeMono (phi.size());     // Size vector
  vector<int>    subMono  (phi.size());     // Subscripts
  vector<double> phiMono  (phi.size());     // Current shape parameters
  vector<vector<double>> PostMono;          // Storage array for the posterior
  vector<vector<long double>> PostMonoExp;  // Storage array for the exponent

  // Declarations for the Multimodal posterior array [Multi]
  // Storage order: PostMulti = [phi1, ..., phiN, pi1, ..., piN]
  // int            nElMulti = 1;               // Number of elements
  // vector<int>    sizeMulti                   // Size vector
  // vector<int>    subMulti                    // Subscripts
  // vector<double> phiMulti  (phi.size());     // Current shape parameters
  // vector<double> PostMulti;                  // Storage array for the posterior
  // vector<long double> PostMultiExp;          // Storage array for the exponent

  Literature     lit(model);

  // Calculate Monomodal posterior elements ------------------------------------
  for (vector<int>::size_type i = 0; i < phi.size(); i++) {
    sizeMono[i] = phi[i].N;
    nElMono *= phi[i].N;
  }

  PostMono.resize(nElMono, vector<double> (data.N));
  PostMonoExp.resize(nElMono, vector<long double> (data.N));

  for (int ind = 0; ind < nElMono; ind++) {
    // From linear index to subscipts
    ind2sub(sizeMono, ind,  subMono);

    // Creating the current combination of values
    for (vector<int>::size_type j = 0; j < phi.size(); j++) {
      phiMono[j] = phi[j].vec[subMono[j]];
    }

    // Cycle over the data
    for (size_t n = 0; n < data.N; n++) {
      // Calculate only the exponent for numerical purpouses
      PostMonoExp[ind][n] = -0.5 * pow( ((data.vt[n]
                          - lit.CalculateVt(data.dv[n], phiMono))
                          / data.sigma[n]), 2 );

    }
  }

  // Rescaling and Exponentiating
  // auto max = max_element(begin(PostMonoExp), end(PostMonoExp));

  for (int ind = 0; ind < nElMono; ind++) {
    for (size_t n = 0; n < data.N; n++) {
      PostMono[ind][n] = exp(PostMonoExp[ind][n]); // - *max);
    }
  }

  // Calculate Multimodal posterior --------------------------------------------





  return PostMono;
}


void Posterior::cumprod(const vector<int>& size,  vector<int>& k) {
  // cumprod Cumulative product

  for (vector<int>::size_type i = 0; i < size.size(); i++) {
    k[i] = 1;
    for (vector<int>::size_type j = 0; j <= i; j++) {
      k[i] *= size[j];
    }
  }
}

void Posterior::ind2sub(const vector<int>& size, int ind,  vector<int>& sub) {
  // IND2SUB Multiple subscripts from linear index.
  //   IND2SUB is used to determine the equivalent subscript values
  //   corresponding to a given single index into an array.
  //
  //   SUB = IND2SUB(SIZE, IND) returns an array of N subscripts
  //   SUB = [S1, S2, ..., SN] containing the subscripts equivalent
  //   to IND for an array of size SIZE.

  vector<int> k (size.size());
  int temp1, temp2;


  cumprod(size, k);

  for (vector<int>::size_type i = k.size(); i > 0; i--) {
    temp1 = ind % k[i-1];
    temp2 = (ind - temp1) / k[i-1];

    sub[i] = temp2;
    ind = temp1;
  }
  sub[0] = ind;
}


void Posterior::sub2ind(const vector<int>& size, const vector<int>& sub,  int& ind) {
  // SUB2IND Linear index from multiple subscripts.
  //   SUB2IND is used to determine the equivalent single index
  //   corresponding to a given set of subscript values.
  //
  //   IND = SUB2IND(SIZE ,SUB) returns the linear index equivalent to
  //   the N subscripts contained in SUB = [S1, S2, ..., SN] for an
  //   array of size SIZE.

  vector<int> k (size.size());
  int         temp;


  cumprod(size,  k);

  temp = size[0];
  for (vector<int>::size_type i = 1; i < k.size(); i++) {
    temp += (sub[i] - 1) * k[i-1];
  }
  ind = temp;
}
