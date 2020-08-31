#include "extendedcauchymatrix.h"
#include <exception>
#include <iomanip>

void ExtendedCauchyMatrix::initialize()
{
  for (unsigned int i = 0; i < getRows(); ++i)
    m_generatorCol.emplace_back(GeneratorElement(i, 1));
  for (unsigned int i = 0; i < getCols(); ++i)
    m_generatorRow.emplace_back(GeneratorElement(i+getRows(), 1));
}

ExtendedCauchyMatrix::ExtendedCauchyMatrix(const unsigned int rows, const unsigned int cols, const Galois &GF):
  m_GF(GF), m_rows(rows), m_cols(cols)
{
  if (getRows() + getCols() > getGF().getMax())
    throw std::runtime_error(
      "No Cauchy matrices of this size in the provided GF"
    );
  initialize();
}

ExtendedCauchyMatrix::ExtendedCauchyMatrix(const ExtendedCauchyMatrix::GeneratorVector &row, unsigned int rows, const Galois &GF):
  m_GF(GF), m_rows(rows), m_cols(static_cast<unsigned int>(row.size())),
  m_generatorRow(row)
{
  std::vector<bool> used(GF.getMax(), false);
  if (getRows() + getCols() > getGF().getMax())
    throw std::runtime_error(
      "No Cauchy matrices of this size in the provided GF"
    );
  for (auto& el: getGeneratorRow()) {
    if (used[el.first]) throw std::runtime_error("Uniquines violated.");
    used[el.first] = true;
  }
  for (unsigned int el = 0; el < GF.getMax() && getGeneratorCol().size() < getRows(); ++el) {
    if (used[el]) continue;
    m_generatorCol.emplace_back(GeneratorElement(el, 1u));
  }
}

ExtendedCauchyMatrix::ExtendedCauchyMatrix(
  const ExtendedCauchyMatrix::GeneratorVector& row,
  const ExtendedCauchyMatrix::GeneratorVector& col, const Galois &GF):
  m_GF(GF), m_rows(static_cast<unsigned int>(col.size())),
  m_cols(static_cast<unsigned int>(row.size())), m_generatorRow(row),
  m_generatorCol(col)
{
  if (getRows() + getCols() > getGF().getMax())
    throw std::runtime_error(
      "No Cauchy matrices of this size in the provided GF"
    );
  std::vector<bool> used(getGF().getMax(), false);
  for (auto& el: getGeneratorRow()) {
    if (used[el.first] == true)
      throw std::runtime_error(
        "Uniquines violated"
      );
    if (el.second == 0u)
      throw std::runtime_error(
        "Zero multiplier"
      );
  }
  for (auto& el: getGeneratorCol()) {
    if (used[el.first] == true)
      throw std::runtime_error(
        "Uniquines violated"
      );
    if (el.second == 0u)
      throw std::runtime_error(
        "Zero multiplier"
      );
  }
}

ExtendedCauchyMatrix::ExtendedCauchyMatrix(const ExtendedCauchyMatrix &other):
  m_GF(other.getGF()), m_rows(other.getRows()), m_cols(other.getCols()),
  m_weight(other.m_weight),
  m_generatorRow(other.getGeneratorRow()),
  m_generatorCol(other.getGeneratorCol())
{
}

ExtendedCauchyMatrix::ExtendedCauchyMatrix(ExtendedCauchyMatrix &&other):
  m_GF(other.getGF()), m_rows(other.getRows()), m_cols(other.getCols()),
  m_weight(other.m_weight)
{
  std::swap(m_generatorRow, other.m_generatorRow);
  std::swap(m_generatorCol, other.m_generatorCol);
}

unsigned int ExtendedCauchyMatrix::getBitmatrixWeight() const
{
  if (m_weight == 0) {
  for (unsigned int row = 0; row < getRows(); ++row)
    for (unsigned int col = 0; col < getCols(); ++col)
      m_weight += getGF().bitmatrixOnes(element(row, col));
  }
  return m_weight;
}

unsigned int ExtendedCauchyMatrix::element(const unsigned int row, const unsigned int col) const
{
  unsigned int y = m_generatorCol[row].first, s = m_generatorCol[row].second;
  unsigned int x = m_generatorRow[col].first, r = m_generatorRow[col].second;
  return m_GF.divide(m_GF.product(r, s), m_GF.sum(x, y));
}

void ExtendedCauchyMatrix::displayGenerator(std::ostream &os) const
{
  os<<std::setfill(' ');
  os<<"Extended Cauchy Matrix ("<<m_GF.getW()<<", "<<getRows()<<", "<<getCols()<<") generators:"<<std::endl;
  os<<"\t x: (";
  for (unsigned int i = 0; i < getCols()-1; ++i)
    os<<std::setw(4)<<m_generatorRow[i].first<<", ";
  os<<std::setw(4)<<m_generatorRow.back().first<<")"<<std::endl;
  os<<"\t r: (";
  for (unsigned int i = 0; i < getCols()-1; ++i)
    os<<std::setw(4)<<m_generatorRow[i].second<<", ";
  os<<std::setw(4)<<m_generatorRow.back().second<<")"<<std::endl;
  os<<"\t y: (";
  for (unsigned int i = 0; i < getRows()-1; ++i)
    os<<std::setw(4)<<m_generatorCol[i].first<<", ";
  os<<std::setw(4)<<m_generatorCol.back().first<<")"<<std::endl;
  os<<"\t s: (";
  for (unsigned int i = 0; i < getRows()-1; ++i)
    os<<std::setw(4)<<m_generatorCol[i].second<<", ";
  os<<std::setw(4)<<m_generatorCol.back().second<<")"<<std::endl;
}

bool ExtendedCauchyMatrix::operator<(const ExtendedCauchyMatrix &other) const
{
  if (other.getRows() != getRows() || other.getCols() != getCols())
    throw std::runtime_error("Two ECRS matrices are incomparable due to dimensions mismatch.");
  if (other.getGF().getW() != getGF().getW())
    throw std::runtime_error("Two ECRS matrices are incomparable due to GF mismatch.");
  return getBitmatrixWeight() < other.getBitmatrixWeight();
}

ExtendedCauchyMatrix &ExtendedCauchyMatrix::operator=(const ExtendedCauchyMatrix &other)
{
  if (other.getRows() != getRows() || other.getCols() != getCols())
    throw std::runtime_error("Can not reassign a ECRS due to size mismatch.");
  if (other.getGF().getW() != getGF().getW())
    throw std::runtime_error("Can not reassign a ECRS due to GF mismatch.");
  m_generatorRow = other.m_generatorRow;
  m_generatorCol = other.m_generatorCol;
  m_weight = other.m_weight;
  return *this;
}

ExtendedCauchyMatrix &ExtendedCauchyMatrix::operator=(ExtendedCauchyMatrix &&other)
{
  if (other.getRows() != getRows() || other.getCols() != getCols())
    throw std::runtime_error("Can not move a ECRS due to size mismatch.");
  if (other.getGF().getW() != getGF().getW())
    throw std::runtime_error("Can not move a ECRS due to GF mismatch.");
  std::swap(m_generatorRow, other.m_generatorRow);
  std::swap(m_generatorCol, other.m_generatorCol);
  std::swap(m_weight, other.m_weight);
  return *this;
}

std::ostream &operator<<(std::ostream &os, const ExtendedCauchyMatrix &matrix)
{
  os<<std::setfill(' ');
  os<<"Extended Cauchy Matrix ("<<matrix.getGF().getW()<<", "
    <<matrix.getRows()<<", "<<matrix.getCols()<<"):"<<std::endl;
  for (unsigned int row = 0; row < matrix.getRows(); ++row) {
    os<<"\t";
    for (unsigned int col = 0; col < matrix.getCols(); ++col) {
      os<<std::setw(4)<<matrix.element(row, col);
    }
    os<<std::endl;
  }
  os<<"Number of ones: "<<matrix.getBitmatrixWeight()<<std::endl;
  matrix.displayGenerator(os);
  return os;
}
