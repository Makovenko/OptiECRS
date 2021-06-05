#include "optimizedirectionalcost.h"

void OptimizeDirectionalCost::initializeCache(const ExtendedCauchyMatrix& m_initial)
{
  const ExtendedCauchyMatrix::GeneratorVector& given(
    m_direction==OptimizeDirection::Direction::Horizontal?
    m_initial.getGeneratorCol():m_initial.getGeneratorRow()
  );
  unsigned int valuesGiven = static_cast<unsigned int>(given.size());

  m_used.resize(m_GF->getMax());
  for (auto& el: given)
    m_used[el.first] = true;

  m_cache.resize(
    m_GF->getMax(),
    std::vector<unsigned int>(m_GF->getMax(), 0)
  );

  std::vector<unsigned int> precalc(valuesGiven, 0);
  for (unsigned int el = 0; el < m_GF->getMax(); ++el) {
    for (unsigned int i = 0; i < valuesGiven; ++i)
      precalc[i] = m_GF->divide(
        given[i].second, m_initial.getGF()->sum(given[i].first, el)
      );
    if (m_used[el]) continue;
    for (unsigned int mul = 1; mul < m_GF->getMax(); ++mul) {
      for (unsigned int i = 0; i < valuesGiven; ++i) {
        m_cache[el][mul] += m_GF->bitmatrixOnes(
          m_GF->product(mul, precalc[i])
        );
      }
    }
  }
}

void OptimizeDirectionalCost::addElement(const ExtendedCauchyMatrix::GeneratorElement& element) {
  m_used[element.first] = true;
  for (unsigned int el = 0; el < m_GF->getMax(); ++el) {
    if (m_used[el]) continue;
    unsigned int precalc = m_GF->divide(element.second, m_GF->sum(element.first, el));
    for (unsigned int mul = 1; mul < m_GF->getMax(); ++mul) {
        m_cache[el][mul] += m_GF->bitmatrixOnes(m_GF->product(mul, precalc));
    }
  }
}

void OptimizeDirectionalCost::removeElement(const ExtendedCauchyMatrix::GeneratorElement &element)
{
  for (unsigned int el = 0; el < m_GF->getMax(); ++el) {
    if (m_used[el]) continue;
    unsigned int precalc = m_GF->divide(element.second, m_GF->sum(element.first, el));
    for (unsigned int mul = 1; mul < m_GF->getMax(); ++mul) {
        m_cache[el][mul] -= m_GF->bitmatrixOnes(m_GF->product(mul, precalc));
    }
  }
  m_used[element.first] = false;
}

OptimizeDirectionalCost::OptimizeDirectionalCost(
  const ExtendedCauchyMatrix& given,
  const OptimizeDirection::Direction &direction
):
  m_GF(given.getGF()),
  m_direction(direction),
  m_valuesNeeded(
    m_direction==OptimizeDirection::Direction::Horizontal?given.getCols():given.getRows()
  )
{
  initializeCache(given);
}

unsigned int OptimizeDirectionalCost::getExtensionCost(const ExtendedCauchyMatrix::GeneratorElement &element)
{
  std::priority_queue<unsigned int> queue;
  for (unsigned int el = 0; el < m_GF->getMax(); ++el) {
    if (m_used[el] || el == element.first) continue;
    unsigned int precalc = m_GF->divide(element.second, m_GF->sum(el, element.first));
    unsigned int minCost = std::numeric_limits<unsigned int>::max();
    for (unsigned int mul = 1; mul < m_GF->getMax(); ++mul) {
      auto cost = m_cache[el][mul] + m_GF->bitmatrixOnes(m_GF->product(mul, precalc));
      minCost = std::min(minCost, cost);
    }
    queue.push(minCost);
    if (queue.size() > m_valuesNeeded) queue.pop();
  }

  unsigned int res = 0;
  while (!queue.empty()) {
    res += queue.top();
    queue.pop();
  }
  return res;
}
