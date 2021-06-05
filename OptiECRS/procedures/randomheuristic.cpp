#include <vector>
#include <random>
#include "randomheuristic.h"
#include "optimizeconsecutive.h"

RandomHeuristic::RandomHeuristic(
    const ExtendedCauchyMatrix &initial,
    unsigned int N
):
  m_initial(initial), m_result(m_initial), m_N(N),
  m_started(std::chrono::high_resolution_clock::now()),
  m_finished(m_started)
{
}

ExtendedCauchyMatrix RandomHeuristic::run(const double timelimit)
{
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<unsigned int> disValue(
    0, m_initial.getGF()->getMax()-1
  );
  std::uniform_int_distribution<unsigned int> disMult(
    1, m_initial.getGF()->getMax()-1
  );
  std::uniform_int_distribution<unsigned int> disPosition(
    0, m_initial.getCols()-1
  );


  auto current = m_initial;
  m_timelimit = timelimit;
  m_started = std::chrono::high_resolution_clock::now();
  current = OptimizeConsecutive(current).run();

  for (unsigned int i = 0; i < m_N; ++i) {
    if (shouldTerminate()) break;
    std::vector<bool> used(m_initial.getGF()->getMax(), false);
    for (auto& e: current.getGeneratorCol()) used[e.first] = true;
    for (auto& e: current.getGeneratorRow()) used[e.first] = true;
    unsigned int pos = disPosition(gen);
    unsigned int value;
    used[current.getGeneratorRow()[pos].first] = false;
    do value = disValue(gen); while (used[value]);
    auto candidate = current;
    candidate.setRowGenerator_unsafe(pos, {value, disMult(gen)});
    candidate = OptimizeConsecutive(candidate).run();
    if (candidate < current) {
      current = candidate;
      std::cerr<<"Found a better candidate on iteration "<<(i+1)<<std::endl;
      std::cerr<<"New best known value is "<<current.getBitmatrixWeight()<<std::endl;
      std::cerr<<"New best known matrix is "<<current<<std::endl;
      std::cerr<<"Timestamp is "<<
        std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::high_resolution_clock::now() - m_started
        ).count() << " seconds."<<std::endl;
    }
  }
  m_result = current;
  m_finished = std::chrono::high_resolution_clock::now();
  return m_result;
}

bool RandomHeuristic::shouldTerminate() const {
  auto dt = std::chrono::high_resolution_clock::now() - m_started;
  return std::chrono::duration_cast<std::chrono::duration<double>>(dt).count() > m_timelimit;
}

double RandomHeuristic::getTimeInMS() const
{
  return std::chrono::duration_cast<std::chrono::duration<double>>(
      m_finished - m_started
    ).count();
}

