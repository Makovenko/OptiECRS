#ifndef RANDOMHEURISTIC_H
#define RANDOMHEURISTIC_H
#include "../matrix/extendedcauchymatrix.h"
#include <chrono>

/*
 * This is a simple random search heuristic. Given a ECRS matrix M and a number
 * of repetitions it procceds as this:
 * 1. Optimize the given ECRS matrix using OptimizeConsecutive procedure. Let
 *    M be the result of the optimization procedure.
 * 2. Randomly select a row generator element of M, replace it with another
 *    random generator element. Let this matrix be M'
 * 3. Optimize M' using OptimizeConsecutive. If M' < M, let M be M'. Go to step
 *    2, unless the limit of repetitions has been reached.
 */

class RandomHeuristic
{
  ExtendedCauchyMatrix m_initial, m_result;
  unsigned int m_N;
  double m_timelimit;

  std::chrono::high_resolution_clock::time_point m_started, m_finished;
  bool shouldTerminate() const;

public:
  RandomHeuristic(const ExtendedCauchyMatrix& m_initial, unsigned int N = 1000);
  ExtendedCauchyMatrix run(const double timelimit = 30);
  double getTimeInMS() const;
};

#endif // RANDOMHEURISTIC_H
