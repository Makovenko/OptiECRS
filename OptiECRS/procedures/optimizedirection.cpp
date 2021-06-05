#include "optimizedirection.h"
#include <queue>
#include <vector>
#include <limits>
#include <numeric>

std::priority_queue<
  std::pair<unsigned int, ExtendedCauchyMatrix::GeneratorElement>
> OptimizeDirection::generate_queue(
  const ExtendedCauchyMatrix &matrix,
  const OptimizeDirection::Direction &direction
) {
  ExtendedCauchyMatrix::GeneratorVector given;
  unsigned int valuesNeeded, valuesGiven;
  switch(direction) {
    case Direction::Horizontal:
      given = matrix.getGeneratorCol();
      valuesGiven = matrix.getRows();
      valuesNeeded = matrix.getCols();
      break;
    case Direction::Vertical:
      given = matrix.getGeneratorRow();
      valuesGiven = matrix.getCols();
      valuesNeeded = matrix.getRows();
      break;
    default:
      valuesNeeded = 0;
      valuesGiven  = 0;
  }

  std::vector<bool> used(matrix.getGF()->getMax(), false);
  for (auto& el: given) used[el.first] = true;
  std::priority_queue<std::pair<unsigned int, ExtendedCauchyMatrix::GeneratorElement>> queue;
  std::vector<unsigned int> precalc(valuesGiven, 0);

  for (unsigned int el = 0; el < matrix.getGF()->getMax(); ++el) {
    if (used[el]) continue;

    for (unsigned int i = 0; i < valuesGiven; ++i)
      precalc[i] = matrix.getGF()->divide(
        given[i].second, matrix.getGF()->sum(given[i].first, el)
      );

    std::pair<unsigned int, unsigned int> best(std::numeric_limits<unsigned int>::max(), 0);
    for (unsigned int mul = 1; mul < matrix.getGF()->getMax(); ++mul) {
      unsigned int cost = 0;
      for (unsigned int i = 0; i < valuesGiven; ++i) {
        cost += matrix.getGF()->bitmatrixOnes(
          matrix.getGF()->product(mul, precalc[i])
        );
      }
      best = min(best, {cost, mul});
    }
    queue.push({best.first, {el, best.second}});
    if (queue.size() > valuesNeeded) queue.pop();
  }
  return queue;
}

ExtendedCauchyMatrix OptimizeDirection::run(const ExtendedCauchyMatrix &matrix, const OptimizeDirection::Direction &direction)
{
  auto queue = generate_queue(matrix, direction);
  ExtendedCauchyMatrix::GeneratorVector result;
  while (!queue.empty()){
    result.emplace_back(queue.top().second);
    queue.pop();
  }

  if (direction == Direction::Horizontal)
    return ExtendedCauchyMatrix(
      result, matrix.getGeneratorCol(), matrix.getGF()
    );
  else
    return ExtendedCauchyMatrix(
      matrix.getGeneratorRow(), result, matrix.getGF()
    );
}

unsigned int OptimizeDirection::run_cost_only(const ExtendedCauchyMatrix &matrix, const OptimizeDirection::Direction &direction)
{
  auto queue = generate_queue(matrix, direction);
  unsigned int result = 0;
  while (!queue.empty()) {
    result += queue.top().first;
    queue.pop();
  }
  return result;
}
