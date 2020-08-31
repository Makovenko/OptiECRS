#include "optimizebranchandbound.h"
#include "optimizedirection.h"
#include "optimizeconsecutive.h"
#include <queue>

void OptimizeBranchAndBound::recursiveCall(const ExtendedCauchyMatrix::GeneratorElement& next)
{
  if (shouldTerminate()) return;
  m_used[next.first] = true;
  m_currentRow.push_back(next);
  m_directionOptimizer.addElement(next);
  if (m_currentRow.size() == m_current.getCols()) {
    m_current = OptimizeDirection::run(
      ExtendedCauchyMatrix(m_currentRow, m_current.getRows(), m_current.getGF()),
      OptimizeDirection::Direction::Vertical
    );
    std::cerr<<"OptimizeBranchAndBound: found a new solution with cost "<<m_current.getBitmatrixWeight()<<std::endl;
    std::cerr<<m_current<<std::endl;
  }
  else {
    unsigned int colsLeft = m_current.getCols() - static_cast<unsigned int>(m_currentRow.size()) - 1;
    using WeightedElement = std::pair<unsigned int, ExtendedCauchyMatrix::GeneratorElement>;
    std::priority_queue<WeightedElement, std::vector<WeightedElement>, std::greater<WeightedElement>> candidates;
    unsigned int maxEl = m_currentRow.size() == 1?2:m_current.getGF().getMax();

//    //for (unsigned int el = next.first + 1; el < maxEl; ++el) {
//    for (unsigned int el = 1; el < maxEl; ++el) {
//      if (m_used[el]) continue;
//      for (unsigned int mul = 1; mul < m_current.getGF().getMax(); ++mul) {
//        auto element = ExtendedCauchyMatrix::GeneratorElement(el, mul);
//        auto cost = m_directionOptimizer.getExtensionCost(element);
//        //if (cost + colsLeft*m_current.getGF().getW()*m_current.getRows() < m_current.getBitmatrixWeight()) {
//        if (cost * m_current.getCols()/static_cast<double>(m_current.getCols()-colsLeft) < m_current.getBitmatrixWeight()) {
//          candidates.push(WeightedElement(cost, element));
//        }
//      }
//    }

    //for (unsigned int el = next.first + 1; el < maxEl; ++el) {
    for (unsigned int el = 1; el < maxEl; ++el) {
      if (m_used[el]) continue;
      auto best = WeightedElement(std::numeric_limits<unsigned int>:: max(), {0, 0});
      for (unsigned int mul = 1; mul < m_current.getGF().getMax(); ++mul) {
        if (shouldTerminate()) return;
        auto element = ExtendedCauchyMatrix::GeneratorElement(el, mul);
        auto cost = m_directionOptimizer.getExtensionCost(element);
        //if (cost + colsLeft*m_current.getGF().getW()*m_current.getRows() < m_current.getBitmatrixWeight())
        if (cost * m_current.getCols()/static_cast<double>(m_current.getCols()-colsLeft) < m_current.getBitmatrixWeight())
          best = std::min(best, {cost, element});
      }
      candidates.push(best);
    }

    while (!candidates.empty()) {
      auto candidate = candidates.top(); candidates.pop();
      //auto bound = candidate.first + colsLeft*m_current.getGF().getW()*m_current.getRows();
      //if (bound >= m_current.getBitmatrixWeight()) break;
      if (candidate.first * m_current.getCols()/static_cast<double>(m_current.getCols()-colsLeft) >= m_current.getBitmatrixWeight()) break;
      recursiveCall(candidate.second);
    }
  }
  m_directionOptimizer.removeElement(next);
  m_currentRow.pop_back();
  m_used[next.first] = false;
}

bool OptimizeBranchAndBound::shouldTerminate() const {
  auto dt = std::chrono::high_resolution_clock::now() - m_started;
  return std::chrono::duration_cast<std::chrono::duration<double>>(dt).count() > m_timelimit;
}

OptimizeBranchAndBound::OptimizeBranchAndBound(const ExtendedCauchyMatrix &initial):
  m_current(initial),
  m_directionOptimizer(ExtendedCauchyMatrix(initial.getRows(), 0, initial.getGF()), OptimizeDirection::Direction::Vertical),
  m_started(std::chrono::high_resolution_clock::now()),
  m_finished(m_started)
{
  m_used.resize(m_current.getGF().getMax(), false);
  m_current = OptimizeConsecutive(m_current).run();
}

ExtendedCauchyMatrix OptimizeBranchAndBound::run(const double timelimit)
{
  m_timelimit = timelimit;
  OptimizeDirectionalCost costCalculator(m_current, OptimizeDirection::Direction::Vertical);
  m_started = std::chrono::high_resolution_clock::now();
  recursiveCall(ExtendedCauchyMatrix::GeneratorElement(0, 1));
  m_finished = std::chrono::high_resolution_clock::now();
  return m_current;
}

double OptimizeBranchAndBound::getTimeInMS() const
{
  return std::chrono::duration_cast<std::chrono::duration<double>>(
    m_finished - m_started
  ).count();
}
