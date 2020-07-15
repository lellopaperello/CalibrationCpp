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

// Brute Force approach
void Posterior::bruteForce (/* arguments */) {

  // Declarations for the Monomodal posterior array [Mono]
  // Storage order: PostMono = [data, phi1, ..., phiN]
  int            nElMono;                   // Number of elements
  vector<int>    sizeMono (phi.size() + 1); // Size vector
  vector<int>    subMono  (phi.size() + 1); // Subscripts
  vector<double> phiMono  (phi.size());     // Current shape parameters
  vector<double> PostMono;                  // Storage array

  Literature     lit(model);

  // Calculate Monomodal posterior elements ------------------------------------
  sizeMono[0] = data.N;
  nElMono = data.N;
  for (vector<int>::size_type i = 0; i < phi.size(); i++) {
    sizeMono[i+1] = phi[i].N;
    nElMono *= phi[i].N;
  }

  PostMono.resize(nElMono);

  for (int ind = 0; ind < nElMono; ind++) {
    // From linear index to subscipts
    ind2sub(sizeMono, ind,  subMono);

    // Creating the current combination of values
    for (vector<int>::size_type j = 0; j < phi.size(); j++) {
      phiMono[j] = phi[j].vec[subMono[j+1]];
    }

    // Calculate only the exponent for numerical purpouses
    PostMono[ind] = (data.vt[subMono[0]]
                  - lit.CalculateVt(data.dv[subMono[0]], phiMono))
                  / data.sigma[subMono[0]];

  }
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
    temp1 = (ind - 1) % (k[i-1]) + 1;
    temp2 = (ind - temp1) / k[i-1] + 1;

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
