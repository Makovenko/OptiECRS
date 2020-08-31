#ifndef OPTIMIZEDIRECTIONALCOST_H
#define OPTIMIZEDIRECTIONALCOST_H
#include "optimizedirection.h"
#include "../matrix/extendedcauchymatrix.h"
#include <vector>

class OptimizeDirectionalCost
{
  std::vector<std::vector<unsigned int>> m_cache;
  std::vector<bool> m_used;
  const Galois m_GF;
  const OptimizeDirection::Direction m_direction;
  const unsigned int m_valuesNeeded;

  void initializeCache(const ExtendedCauchyMatrix& given);
  void updateCache(const ExtendedCauchyMatrix::GeneratorElement &element);

public:
  OptimizeDirectionalCost(const ExtendedCauchyMatrix& given, const OptimizeDirection::Direction& direction);
  unsigned int getExtensionCost(const ExtendedCauchyMatrix::GeneratorElement& el);
  void addElement(const ExtendedCauchyMatrix::GeneratorElement& element);
  void removeElement(const ExtendedCauchyMatrix::GeneratorElement& element);
};

#endif // OPTIMIZEDIRECTIONALCOST_H
