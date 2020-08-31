#include "optimizeconsecutive.h"
#include "optimizedirection.h"

OptimizeConsecutive::OptimizeConsecutive(const ExtendedCauchyMatrix &initial):
  m_initial(initial), m_result(initial)
{
  m_started = std::chrono::high_resolution_clock::now();
  m_finished = m_started;
}

ExtendedCauchyMatrix OptimizeConsecutive::run()
{
  m_started = std::chrono::high_resolution_clock::now();
  // First run: initial direction is horizontal
  for (OptimizeDirection::Direction direction:
    {OptimizeDirection::Direction::Horizontal,
     OptimizeDirection::Direction::Vertical}
  ) {
    auto current = m_initial, best = m_initial;
    do {
      current = best;
      auto better = OptimizeDirection::run(current, direction);
      best = better;
      direction = direction==OptimizeDirection::Direction::Horizontal?
        OptimizeDirection::Direction::Vertical:
        OptimizeDirection::Direction::Horizontal;
    } while (best < current);
    m_result = std::min(m_result, best);
  }
  m_finished = std::chrono::high_resolution_clock::now();
  return m_result;
}

double OptimizeConsecutive::getTimeInMS() const
{
  return std::chrono::duration_cast<std::chrono::duration<double>>(
      m_finished - m_started
    ).count();
}
