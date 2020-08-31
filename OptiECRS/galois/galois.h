#ifndef GALOIS_H
#define GALOIS_H
#include <algorithm>
#include <vector>

/*
 * This is a wrapper class for gf_complete and jerasure libraries routines.
 * All interactions with GF(2^w) in the project should be implemented through
 * this class.
 */

class Galois
{
private:
  unsigned int m_w;
  std::vector<std::vector<unsigned int>> m_sumCache;
  std::vector<std::vector<unsigned int>> m_productCache;
  std::vector<std::vector<unsigned int>> m_divideCache;
  std::vector<unsigned int> m_numberOfOnesCache;

  void precomputeSums();
  void precomputeProducts();
  void precomputeInverse();
  void precomputeDivision();
  void precomputeNumberOfOnes();

public:
  Galois(const unsigned int w);

  inline unsigned int getW() const { return m_w; }
  inline unsigned int getMax() const { return ((unsigned int)1) << m_w; }

  inline unsigned int sum(const unsigned int rhs, const unsigned int lhs) const {
    return m_sumCache[rhs][lhs];
  }
  inline unsigned int product(const unsigned int rhs, const unsigned int lhs) const {
    return m_productCache[rhs][lhs];
  }
  inline unsigned int divide(const unsigned int rhs, const unsigned int lhs) const {
    return m_divideCache[rhs][lhs];
  }
  inline unsigned int invert(const unsigned int rhs) const {
    return m_divideCache[1][rhs];
  }
  inline unsigned int bitmatrixOnes(const unsigned int rhs) const {
    return m_numberOfOnesCache[rhs];
  }

  #ifdef DEBUG_GALOIS
  inline unsigned int sum_uncached(const unsigned int rhs, const unsigned int lhs) const {
    unsigned int res = rhs;
    unsigned int l   = lhs;
    galois_region_xor(reinterpret_cast<char*>(&l), reinterpret_cast<char*>(&res), sizeof(res));
    return res;
  }
  inline unsigned int product_uncached(const unsigned int rhs, const unsigned int lhs) const {
    return static_cast<unsigned int>(galois_single_multiply(static_cast<int>(rhs), static_cast<int>(lhs), static_cast<int>(m_w)));
  }
  inline unsigned int invert_uncached(const unsigned int rhs) const {
    return static_cast<unsigned int>(galois_inverse(static_cast<int>(rhs), static_cast<int>(m_w)));
  }
  inline unsigned int divide_uncached(const unsigned int rhs, const unsigned int lhs) const {
    return static_cast<unsigned int>(galois_single_divide(static_cast<int>(rhs), static_cast<int>(lhs), static_cast<int>(m_w)));
  }

  void comparePerformance() const;
  #endif
};

#endif // GALOIS_H
