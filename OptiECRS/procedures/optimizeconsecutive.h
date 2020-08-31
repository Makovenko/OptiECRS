#ifndef OPTIMIZECONSECUTIVE_H
#define OPTIMIZECONSECUTIVE_H
#include <chrono>
#include "../matrix/extendedcauchymatrix.h"

/*
 * This is a simple consequtive optimization algorithm based on the
 * OptimizeDirectionProcedure. It proceeds by optimizing a given ECRS
 * over alterniting directions. Since there is a choice of an initial
 * optimization direction, it performs two runs: one starting with
 * horizontal direction, another one -- with vertical.
 */

class OptimizeConsecutive
{
  ExtendedCauchyMatrix m_initial;
  ExtendedCauchyMatrix m_result;
  std::chrono::high_resolution_clock::time_point m_started, m_finished;

public:
  OptimizeConsecutive(const ExtendedCauchyMatrix& initial);
  ExtendedCauchyMatrix run();
  double getTimeInMS() const;
};

#endif // OPTIMIZECONSECUTIVE_H
