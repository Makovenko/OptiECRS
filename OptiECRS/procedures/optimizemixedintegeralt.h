#ifndef OPTIMIZEMIXEDINTEGERALT_H
#define OPTIMIZEMIXEDINTEGERALT_H
#include "../matrix/extendedcauchymatrix.h"
#include <vector>
#include <memory>
#include "gurobi_c++.h"

class OptimizeMixedIntegerAlt
{
  const std::shared_ptr<Galois> m_GF;
  GRBEnv m_grbEnv;
  GRBModel m_grbModel;
  std::vector<GRBVar> m_X, m_Y, m_C, m_D, m_U;
  unsigned int m_cols, m_rows;
  
  inline unsigned int tripletToIndexX(const unsigned int first, const unsigned int second, const unsigned int third) {
    return (m_GF->getMax()*first + second)*m_cols + third; }
  inline unsigned int tripletToIndexY(const unsigned int first, const unsigned int second, const unsigned int third) {
  return (m_GF->getMax()*first + second)*m_rows + third; }
  inline unsigned int tripletToIndexCUD(const unsigned int first, const unsigned int second, const unsigned int third) {
  return (m_cols*first + second)*m_rows + third; }

  void initializeVariables(const ExtendedCauchyMatrix& given);
  void initializeConstraints(const ExtendedCauchyMatrix& given);

public:
  OptimizeMixedIntegerAlt(const ExtendedCauchyMatrix& given);
  ExtendedCauchyMatrix run(const double timelimit = 30);
};

#endif // OPTIMIZEMIXEDINTEGERALT_H
