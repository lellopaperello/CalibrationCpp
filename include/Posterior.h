// Header for the Posterior Class.

#ifndef POSTERIOR_H
#define POSTERIOR_H

#include <string>
#include <vector>

Class Posterior
{
  public:

  // Methods
  void cumprod(const vector<int>& size,  vector<int>& k);
  void ind2sub(const vector<int>& size, int ind,  vector<int>& sub);
  void sub2ind(const vector<int>& size, const vector<int>& sub,  int& ind);

};
#endif // POSTERIOR_H
