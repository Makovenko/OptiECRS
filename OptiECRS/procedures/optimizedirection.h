#ifndef OPTIMIZEDIRECTION_H
#define OPTIMIZEDIRECTION_H
#include <queue>
#include "../matrix/extendedcauchymatrix.h"

class OptimizeDirection
{
public:
  enum class Direction{Horizontal, Vertical};
private:
  static std::priority_queue<
    std::pair<unsigned int, ExtendedCauchyMatrix::GeneratorElement>
  > generate_queue(
    const ExtendedCauchyMatrix &matrix,
    const Direction &direction
  );
public:
  static ExtendedCauchyMatrix run(const ExtendedCauchyMatrix& matrix, const Direction& direction);
  static unsigned int run_cost_only(const ExtendedCauchyMatrix& matrix, const Direction& direction);
};

#endif // OPTIMIZEDIRECTION_H
