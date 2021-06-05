#ifndef EXTENDEDCAUCHYMATRIX_H
#define EXTENDEDCAUCHYMATRIX_H
#include "../galois/galois.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

class ExtendedCauchyMatrix
{
public:
  using GeneratorElement = std::pair<unsigned int, unsigned int>;
  using GeneratorVector  = std::vector<GeneratorElement>;

private:
  const std::shared_ptr<Galois> m_GF;
  const unsigned int m_rows, m_cols;
  mutable unsigned int m_weight = 0;

  GeneratorVector m_generatorRow;
  GeneratorVector m_generatorCol;

  /* Initializes a standard CRS code as described in the original paper
   * by Blomer et al, 1995
   */
  void initialize();

public:
  ExtendedCauchyMatrix(const unsigned int rows, const unsigned int cols, const std::shared_ptr<Galois> GF);
  ExtendedCauchyMatrix(const GeneratorVector& row, const GeneratorVector& col, const std::shared_ptr<Galois> GF);
  ExtendedCauchyMatrix(const GeneratorVector& row, unsigned int rows, const std::shared_ptr<Galois> GF);
  ExtendedCauchyMatrix(const ExtendedCauchyMatrix& other);
  ExtendedCauchyMatrix(ExtendedCauchyMatrix&& other);

  inline unsigned int getRows() const { return m_rows; }
  inline unsigned int getCols() const { return m_cols; }
  inline const std::shared_ptr<Galois> getGF()  const { return m_GF; }

  unsigned int getBitmatrixWeight() const;

  GeneratorVector getGeneratorRow() const { return m_generatorRow; }
  GeneratorVector getGeneratorCol() const { return m_generatorCol; }
  inline unsigned int element(const unsigned int row, const unsigned int col) const;

  void displayGenerator(std::ostream& os) const;

  bool operator<(const ExtendedCauchyMatrix& other) const;
  ExtendedCauchyMatrix& operator=(const ExtendedCauchyMatrix& other);
  ExtendedCauchyMatrix& operator=(ExtendedCauchyMatrix&& other);

  void setRowGenerator_unsafe(unsigned int position, const GeneratorElement& el) {
    m_generatorRow[position] = el;
    m_weight = 0;
  }
};

std::ostream& operator<<(std::ostream& os, const ExtendedCauchyMatrix& matrix);

#endif // EXTENDEDCAUCHYMATRIX_H
