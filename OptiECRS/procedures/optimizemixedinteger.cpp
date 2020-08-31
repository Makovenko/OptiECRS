#include "optimizemixedinteger.h"
#include <iostream>

OptimizeMixedInteger::OptimizeMixedInteger(const ExtendedCauchyMatrix& given):
  m_GF(given.getGF()), m_grbEnv(), m_grbModel(m_grbEnv) {
    initializeVariables(given);
    initializeConstraints(given);
}

void OptimizeMixedInteger::initializeVariables(const ExtendedCauchyMatrix& given) {
  for (unsigned int i = 0; i < m_GF.getMax()*m_GF.getMax(); ++i) {
    m_X.push_back(m_grbModel.addVar(0.0, 1.0, 0.0, GRB_BINARY));
    m_X.back().set(GRB_DoubleAttr_Start, 0.0);
    m_Y.push_back(m_grbModel.addVar(0.0, 1.0, 0.0, GRB_BINARY));
    m_Y.back().set(GRB_DoubleAttr_Start, 0.0);
    for (unsigned int j = 0; j < m_GF.getMax()*m_GF.getMax(); ++j) {
      unsigned int u = indexToFirst(i), v = indexToSecond(i);
      unsigned int p = indexToFirst(j), q = indexToSecond(j);
      unsigned int c = m_GF.bitmatrixOnes(m_GF.divide(m_GF.product(v, q), m_GF.sum(u, p)));
      m_T.push_back(m_grbModel.addVar(0.0, 1.0, c, GRB_CONTINUOUS));
    }
  }
  for (auto& el: given.getGeneratorRow()) {
    unsigned int index = pairToIndex(el.first, el.second);
    m_X[index].set(GRB_DoubleAttr_Start, 1.0);
  }
  for (auto& el: given.getGeneratorCol()) {
    unsigned int index = pairToIndex(el.first, el.second);
    m_Y[index].set(GRB_DoubleAttr_Start, 1.0);
  }
}

void OptimizeMixedInteger::initializeConstraints(const ExtendedCauchyMatrix &given) {
  GRBLinExpr sumOfX, sumOfY;
  for (unsigned int first = 0; first < m_GF.getMax(); ++first) {
    GRBLinExpr sumOfXPart, sumOfYPart;
    for (unsigned int second = 1; second < m_GF.getMax(); ++second) {
      sumOfXPart += m_X[pairToIndex(first, second)];
      sumOfYPart += m_Y[pairToIndex(first, second)];
    }
    m_grbModel.addConstr(sumOfXPart + sumOfYPart, GRB_LESS_EQUAL, 1.0);
    sumOfX += sumOfXPart;
    sumOfY += sumOfYPart;
  }
  m_grbModel.addConstr(sumOfX, GRB_EQUAL, given.getCols());
  m_grbModel.addConstr(sumOfY, GRB_EQUAL, given.getRows());
  
  for (unsigned int i = 0; i < m_GF.getMax()*m_GF.getMax(); ++i) {
    for (unsigned int j = 0; j < m_GF.getMax()*m_GF.getMax(); ++j) {
      m_grbModel.addConstr(m_T[i*m_GF.getMax()*m_GF.getMax() + j], GRB_GREATER_EQUAL, m_X[i] + m_Y[j] - 1);
    }
  }
}

ExtendedCauchyMatrix OptimizeMixedInteger::run(const double timelimit) {
  m_grbModel.set(GRB_DoubleParam_TimeLimit, timelimit);
  m_grbModel.optimize();
  ExtendedCauchyMatrix::GeneratorVector row, col;
  for (unsigned int i = 0; i < m_GF.getMax()*m_GF.getMax(); ++i) {
    ExtendedCauchyMatrix::GeneratorElement el = {indexToFirst(i), indexToSecond(i)};
    if (m_X[i].get(GRB_DoubleAttr_X) > 0.5) {
      row.push_back(el);
    }
    if (m_Y[i].get(GRB_DoubleAttr_X) > 0.5) {
      col.push_back(el);
    }
  }
  return ExtendedCauchyMatrix(row, col, m_GF);
}
