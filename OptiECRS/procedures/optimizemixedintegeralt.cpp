#include "optimizemixedintegeralt.h"
#include <iostream>

OptimizeMixedIntegerAlt::OptimizeMixedIntegerAlt(const ExtendedCauchyMatrix& given):
m_GF(given.getGF()), m_grbEnv(), m_grbModel(m_grbEnv) {
    m_cols = given.getCols();
    m_rows = given.getRows();
    initializeVariables(given);
    initializeConstraints(given);
}

void OptimizeMixedIntegerAlt::initializeVariables(const ExtendedCauchyMatrix& given) {
  for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
    for (unsigned int v = 0; v < m_GF->getMax(); ++v) {
      for (unsigned int i = 0; i < given.getCols(); ++i) {
        m_X.push_back(m_grbModel.addVar(0.0, 1.0, 0.0, GRB_BINARY));
      }
    }
  }
  for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
    for (unsigned int v = 0; v < m_GF->getMax(); ++v) {
      for (unsigned int i = 0; i < given.getRows(); ++i) {
        m_Y.push_back(m_grbModel.addVar(0.0, 1.0, 0.0, GRB_BINARY));
      }
    }
  }

  for (unsigned int i = 0; i < given.getCols(); ++i) {
    auto el = given.getGeneratorRow()[i];
    m_X[tripletToIndexX(el.first, el.second, i)].set(GRB_DoubleAttr_Start, 1.0);
  }
  for (unsigned int j = 0; j < given.getRows(); ++j) {
    auto el = given.getGeneratorCol()[j];
    m_Y[tripletToIndexY(el.first, el.second, j)].set(GRB_DoubleAttr_Start, 1.0);
  }
  
  for (unsigned int k = 0; k < m_GF->getMax()*given.getCols()*given.getRows(); ++k) {
    m_C.push_back(m_grbModel.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS));
    m_D.push_back(m_grbModel.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS));
  }
  for (unsigned int t = 0; t < m_GF->getMax(); ++t) {
    for (unsigned int k = 0; k < given.getCols()*given.getRows(); ++k) {
      m_U.push_back(m_grbModel.addVar(
        0.0, 1.0,
        t==0?given.getCols()*given.getRows()*m_GF->getW()*m_GF->getW():m_GF->bitmatrixOnes(t),
        GRB_CONTINUOUS
      ));
    }
  }
}

void OptimizeMixedIntegerAlt::initializeConstraints(const ExtendedCauchyMatrix &given) {
  for (unsigned int k = 0; k < given.getCols(); ++k) {
    GRBLinExpr xSum;
    for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
      for (unsigned int v = 1; v < m_GF->getMax(); ++v) {
        xSum += m_X[tripletToIndexX(u, v, k)];
      }
    }
    m_grbModel.addConstr(xSum, GRB_EQUAL, 1.0);
  }

  for (unsigned int k = 0; k < given.getRows(); ++k) {
    GRBLinExpr ySum;
    for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
      for (unsigned int v = 1; v < m_GF->getMax(); ++v) {
        ySum += m_Y[tripletToIndexY(u, v, k)];
      }
    }
    m_grbModel.addConstr(ySum, GRB_EQUAL, 1.0);
  }
  
  for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
    GRBLinExpr sum;
    for (unsigned int v = 1; v < m_GF->getMax(); ++v) {
      for (unsigned int k = 0; k < given.getCols(); ++k) {
        sum += m_X[tripletToIndexX(u, v, k)];
      }
      for (unsigned int k = 0; k < given.getRows(); ++k) {
        sum += m_Y[tripletToIndexY(u, v, k)];
      }
    }
    m_grbModel.addConstr(sum, GRB_LESS_EQUAL, 1.0);
  }

  for (unsigned int i = 0; i < given.getCols(); ++i) {
    for (unsigned int j = 0; j < given.getRows(); ++j) {
      for (unsigned int v = 1; v < m_GF->getMax(); ++v) {
        for (unsigned int q = 1; q < m_GF->getMax(); ++q) {
          unsigned int t = m_GF->product(v, q);
          GRBLinExpr sum;
          for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
            sum += m_X[tripletToIndexX(u, v, i)];
            sum += m_Y[tripletToIndexY(u, q, j)];
          }
          m_grbModel.addConstr(m_C[tripletToIndexCUD(t, i, j)], GRB_GREATER_EQUAL, sum - 1);
        }
      }
    }
  }

  for (unsigned int i = 0; i < given.getCols(); ++i) {
    for (unsigned int j = 0; j < given.getRows(); ++j) {
      for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
        for (unsigned int p = 0; p < m_GF->getMax(); ++p) {
          unsigned int t = m_GF->sum(u, p);
          GRBLinExpr sum;
          for (unsigned int v = 1; v < m_GF->getMax(); ++v) {
            sum += m_X[tripletToIndexX(u, v, i)];
            sum += m_Y[tripletToIndexY(p, v, j)];
          }
          m_grbModel.addConstr(m_D[tripletToIndexCUD(t, i, j)], GRB_GREATER_EQUAL, sum - 1);
        }
      }
    }
  }

  for (unsigned int i = 0; i < given.getCols(); ++i) {
    for (unsigned int j = 0; j < given.getRows(); ++j) {
      for (unsigned int t = 1; t < m_GF->getMax(); ++t) {
        for (unsigned int p = 1; p < m_GF->getMax(); ++p) {
          unsigned int u = m_GF->product(t, p);
          auto U = m_U[tripletToIndexCUD(t, i, j)];
          auto C = m_C[tripletToIndexCUD(u, i, j)];
          auto D = m_D[tripletToIndexCUD(p, i, j)];
          m_grbModel.addConstr(U, GRB_GREATER_EQUAL, C + D - 1);
        }
      }
    }
  }
}

ExtendedCauchyMatrix OptimizeMixedIntegerAlt::run(const double timelimit) {
  m_grbModel.set(GRB_DoubleParam_TimeLimit, timelimit);
  m_grbModel.optimize();
  ExtendedCauchyMatrix::GeneratorVector row, col;
  for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
    for (unsigned int v = 0; v < m_GF->getMax(); ++v) {
      for (unsigned int i = 0; i < m_cols; ++i) {
        if (m_X[tripletToIndexX(u, v, i)].get(GRB_DoubleAttr_X) > 0.5) {
          row.push_back({u, v});
        }
      }
    }
  }
  for (unsigned int u = 0; u < m_GF->getMax(); ++u) {
    for (unsigned int v = 0; v < m_GF->getMax(); ++v) {
      for (unsigned int i = 0; i < m_rows; ++i) {
        if (m_Y[tripletToIndexY(u, v, i)].get(GRB_DoubleAttr_X) > 0.5) {
          col.push_back({u, v});
        }
      }
    }
  }

  return ExtendedCauchyMatrix(row, col, m_GF);
}
