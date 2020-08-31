#ifndef OPTIMIZEMIXEDINTEGERALT_H
#define OPTIMIZEMIXEDINTEGERALT_H
#include "../matrix/extendedcauchymatrix.h"
#include <vector>
#include <memory>
#include "gurobi_c++.h"

class OptimizeMixedIntegerAlt
{
  const Galois m_GF;
  GRBEnv m_grbEnv;
  GRBModel m_grbModel;
  std::vector<GRBVar> m_X, m_Y, m_C, m_D, m_U;
  unsigned int m_cols, m_rows;
  
  inline unsigned int indexToFirst(const unsigned int index) const { return index % m_GF.getMax(); }
  inline unsigned int indexToSecond(const unsigned int index) const { return index / m_GF.getMax(); }
  inline unsigned int pairToIndex(const unsigned int first, const unsigned int second) { return second*m_GF.getMax() + first; }

  void initializeVariables(const ExtendedCauchyMatrix& given);
  void initializeConstraints(const ExtendedCauchyMatrix& given);

public:
  OptimizeMixedIntegerAlt(const ExtendedCauchyMatrix& given);
  ExtendedCauchyMatrix run(const double timelimit = 30);
};

#endif // OPTIMIZEMIXEDINTEGERALT_H
