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

// Constructor
Posterior::Posterior(const data_t& DATA, const vector<phi_t>& PHI,
                     const vector<double>& PI, const string& MODEL) {
  data = DATA;
  phi = PHI;
  pi = PI;
  model = MODEL;
};

// Brute Force approach
vector<double> Posterior::bruteForce () {

  // General Declarations
  Literature lit(model);
  const long double realmin = numeric_limits<long double>::min();

  // Declarations for the Monomodal posterior array [Mono]
  // Storage order: PostMono = [phi1, ..., phiN | nData]
  int                         nElMono = 1;            // Number of elements
  vector<int>                 sizeMono (phi.size());  // Size vector
  vector<int>                 subMono  (phi.size());  // Subscripts
  vector<double>              phiMono  (phi.size());  // Current shape parameters

  long double                 maxMonoExp = -1e10;     // Maximum of the posterior
  long double                 minMonoExp = 1e10;      // Minimum of the posterior

  // Calculate Monomodal posterior elements ------------------------------------
  for (vector<int>::size_type i = 0; i < phi.size(); i++) {
    sizeMono[i] = phi[i].vec.size();
    nElMono *= phi[i].vec.size();
  }

  // Storage array for the Monomodal Posterior
  vector<vector<long double>> PostMonoExp (nElMono, vector<long double> (data.N));
  vector<vector<long double>> PostMono    (nElMono, vector<long double> (data.N));

  // Cycle over all the possible combination
  for (int ind = 0; ind < nElMono; ind++) {
    // From linear index to subscipts
    ind2sub(sizeMono, ind,  subMono);

    // Creating the current combination of values
    for (vector<int>::size_type j = 0; j < phi.size(); j++) {
      phiMono[j] = phi[j].vec[subMono[j]];
    }

    // Cycle over the data
    for (int n = 0; n < data.N; n++) {
      // Calculate only the exponent for numerical purpouses
      PostMonoExp[ind][n] = -0.5 * pow( ((data.vt[n]
                          - lit.CalculateVt(data.dv[n], phiMono))
                          / data.sigma[n]), 2 );

      // Updating maximum and minimum
      if (PostMonoExp[ind][n] > maxMonoExp) {
        maxMonoExp = PostMonoExp[ind][n];
      }else if (PostMonoExp[ind][n] < minMonoExp) {
        minMonoExp = PostMonoExp[ind][n];
      }
    }
  }

  /* -------------------------------------------------------------------------- */
  // Checking the number of modes (Ks)
  int K = 1;
  for (vector<int>::size_type i = 0; i < phi.size(); i++) {
    K *= phi[i].K;
  }

  if (K == 1) { // Monomodal posterior
    vector<long double> PostExp   (nElMono, 0);
    vector<double>      Posterior (nElMono, 0);
    long double maxExp = -1e10;

    for (int ind = 0; ind < nElMono; ind++) {
      // Summation over the data
      for (int n = 0; n < data.N; n++) {
        PostExp[ind] += PostMonoExp[ind][n];
      }
      if (PostExp[ind] > maxExp) {
        maxExp = PostExp[ind];
      }
    }

    // Rescaling and Exponentiating
    for (int ind = 0; ind < nElMono; ind++) {
      Posterior[ind] = exp(PostExp[ind] - maxExp);
    }

    return Posterior;
  } else { // Multimodal posterior
    // Rescaling and Exponentiating
    long double base = exp(log(realmin) / (abs(maxMonoExp - minMonoExp)));
    for (int ind = 0; ind < nElMono; ind++) {
      for (int n = 0; n < data.N; n++) {
        PostMono[ind][n] = exp(PostMonoExp[ind][n] - maxMonoExp);
      }
    }

    // Declarations for the Multimodal posterior array [Multi]
    // Storage order: PostMulti = [phi1_1, ..., phi1_k1, ...,
    //                             phiN_1, ..., phiN_kN, ...,
    //                             pi1, ..., piK]
    // With K = k1 * ... * kN.

    int                 nElMulti = 1;           // Number of elements
    vector<int>         sizeMulti;              // Global size vector
    vector<int>         sizeK (phi.size());     // K size vector

    long double         maxMultiExp = -1e10;    // Maximum of the posterior
    long double         minMultiExp = 1e10;     // Minimum of the posterior

    // Calculate Multimodal posterior --------------------------------------------
    for (vector<int>::size_type i = 0; i < phi.size(); i++) {
      for (int k = 0; k < phi[i].K; k++) {
        // Building the size container of the PHIs.
        sizeMulti.push_back(phi[i].vec.size());

        // Total number of combinations
        nElMulti *= phi[i].vec.size();
      }
      // Building the size container of each PHI's K.
      sizeK[i] = phi[i].K;
    }

    // Building the size container of the global indipendent Ks.
    vector<int> sizeGlobalK (K-1, pi.size());
    sizeMulti.insert(sizeMulti.end(), sizeGlobalK.begin(), sizeGlobalK.end());
    vector<int> subMulti (sizeMulti.size());

    // Total number of combinations
    nElMulti *= ((K-1) * pi.size());

    // Storage array for the Multimodal Posterior
    vector<double>      PostMulti    (nElMulti, 0);
    vector<long double> PostMultiExp (nElMulti, 0);

    // Cycle over all the possible combination
    for (int ind = 0; ind < nElMulti; ind++) {
      // From linear index to subscipts
      ind2sub(sizeMulti, ind,  subMulti);

      // Extract coefficeints and enforce normalization constraint
      vector<double> curr_pi;
      for (vector<int>::size_type i = subMulti.size()-K+1; i < subMulti.size(); i++) {
        curr_pi.push_back(pi[subMulti[i]]);
      }
      vector<double> coeff = normConstraint(curr_pi);

      // Generate container for the current parameter (PHI) indices
      vector<vector<int>> phiInd (phi.size());
      for (size_t i = 0; i < phi.size(); i++) {
        int phiInd_i = 0;
        for (size_t j = 1; j <= i; j++) {
          phiInd_i += sizeK[j-1];
        }

        phiInd[i] = vector<int> (subMulti.begin() + phiInd_i,
                                 subMulti.begin() + phiInd_i + sizeK[i]);
      }

      // Inner Summation (k)
      vector<long double> innerSum (data.N, 0);
      for (int k = 0; k < K; k++) {
        // Selecting the Posterior element index corresponding to k
        vector<int> subK (phi.size());
        ind2sub(sizeK, k,  subK);

        vector<int> phiSub (phi.size());
        for (size_t i = 0; i < phi.size(); i++) {
          phiSub[i] = phiInd[i][subK[i]];
        }
        int innerInd;
        sub2ind(sizeMono, phiSub,  innerInd);

        // Performing the summation;
        vector<long double> rhs = (PostMono[innerInd] * coeff[k]);
        innerSum = innerSum + rhs;
      }

      // Outer Summation (n) ~ Using the logarithm for numerical purpouses.
      for (int n = 0; n < data.N; n++) {
        PostMultiExp[ind] += (log(innerSum[n])); // / log(base));
      }
      // Updating maximum and minimum
      if (PostMultiExp[ind] > maxMultiExp) {
        maxMultiExp = PostMultiExp[ind];
      }else if (PostMultiExp[ind] < minMultiExp) {
        minMultiExp = PostMultiExp[ind];
      }
    }

    for (int ind = 0; ind < nElMulti; ind++) {
      PostMulti[ind] = exp(PostMultiExp[ind] - maxMultiExp);
    }

    return PostMulti;
  }
}


vector<int> Posterior::cumprod(const vector<int>& size) {
  // cumprod Cumulative product

  vector<int> k (size.size());
  for (vector<int>::size_type i = 0; i < size.size(); i++) {
    k[i] = 1;
    for (vector<int>::size_type j = 0; j <= i; j++) {
      k[i] *= size[j];
    }
  }
  return k;
}

void Posterior::ind2sub(const vector<int>& size, int ind,  vector<int>& sub) {
  // IND2SUB Multiple subscripts from linear index.
  //   IND2SUB is used to determine the equivalent subscript values
  //   corresponding to a given single index into an array.
  //
  //   SUB = IND2SUB(SIZE, IND) returns an array of N subscripts
  //   SUB = [S1, S2, ..., SN] containing the subscripts equivalent
  //   to IND for an array of size SIZE.

  if (size.size() != sub.size()) {
    throw ; // QUALCOSA
  }

  vector<int> k = cumprod(size);

  for (vector<int>::size_type i = sub.size()-1; i > 0; i--) {
    int temp1 = ind % k[i-1];
    int temp2 = (ind - temp1) / k[i-1];

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

  vector<int> k = cumprod(size);
  int temp = sub[0];

  for (vector<int>::size_type i = 0; i < k.size(); i++) {
    temp += sub[i] * k[i-1];
  }
  ind = temp;
}

vector<double> Posterior::normConstraint(vector<double> pi) {
  // Enforcing the normalization constraint on the K-1 indipendent choices of
  // the mixing coefficient to get the K values for which:
  //                  sum(pi_star) = 1

  pi.push_back(1.0);
  vector<double> pi_star (pi.size());

  for (size_t k = 0; k < pi.size(); k++) {
    pi_star[k] = 1;
    for (size_t i = 0; i <= k; i++) {
      pi_star[k] *=  (i == k) ? pi[i] : (1 - pi[i]);
    }
  }
  return pi_star;
}

long double Posterior::findMax(const vector<vector<long double>>& vec, int s1, int s2) {
  vector<long double> maxVec (s2);
  for (int ind = 0; ind < s1; ind++) {
    maxVec[ind] = *max_element(begin(vec[ind]), end(vec[ind]));
  }
  auto maxMono = max_element(begin(maxVec), end(maxVec));

  return *maxMono;
}


//vector<vector<int>> Posterior::combVec (const vector<vector<int>>& X) {
//  // COMBVEC Create all combinations of vectors.
//  //
//  // COMBVEC takes as input a vector X containing any number of vectors, where
//  // each Xi has Ni elements and return a vector of (N1*N2*...) vectors, where
//  // the inner vectors consist of all combinations found by combining one
//  // element from each Xi.
//
//  int                 nComb = 1;        // Number of total combinations
//  int                 nVec = X.size();  // Number of vectors
//  vector<int>         size (X.size());  // Size vector
//  vector<int>         sub (X.size());   // Subscripts vector
//  vector<vector<int>> Y;                // Result
//
//  for (size_t i = 0; i < nVec; i++) {
//    nComb *= X[i].size();
//    size[i] = X[i].size();
//  }
//
//  Y.resize(nComb, vector<int> (nVec));
//  for (size_t ind = 0; ind < nComb; ind++) {
//    ind2sub(size, ind, sub);
//
//    for (size_t i = 0; i < nVec; i++) {
//      Y[ind][i] = X[i][sub[i]];
//    }
//  }
//
//  return Y;
//}
