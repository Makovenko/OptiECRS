#ifndef OPTIMIZEBRANCHANDBOUND_H
#define OPTIMIZEBRANCHANDBOUND_H
#include <chrono>
#include "../matrix/extendedcauchymatrix.h"
#include "optimizedirectionalcost.h"

class OptimizeBranchAndBound
{
  ExtendedCauchyMatrix m_current;
  ExtendedCauchyMatrix::GeneratorVector m_currentRow;
  OptimizeDirectionalCost m_directionOptimizer;
  std::vector<bool> m_used;
  double m_timelimit;

  std::chrono::high_resolution_clock::time_point m_started, m_finished;

  void recursiveCall(const ExtendedCauchyMatrix::GeneratorElement& next);
  bool shouldTerminate() const;

public:
  OptimizeBranchAndBound(const ExtendedCauchyMatrix& initial);
  ExtendedCauchyMatrix run(const double timelimit = 30);
  double getTimeInMS() const;
};

#endif // OPTIMIZEBRANCHANDBOUND_H
