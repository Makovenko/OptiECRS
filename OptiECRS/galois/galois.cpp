#include "galois.h"
#include "../external/bundle/galois.h"
#include "../external/bundle/jerasure.h"
#include <exception>

#ifdef DEBUG_GALOIS
#include <iostream>
#include <iomanip>
#include <chrono>
#endif

void Galois::precomputeSums()
{
  m_sumCache.resize(getMax(), std::vector<unsigned int>(getMax(), 0));
  for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
    for (unsigned int lhs = rhs; lhs < getMax(); ++lhs) {
      auto tar = reinterpret_cast<char*>(&m_sumCache[rhs][lhs]);
      galois_region_xor(reinterpret_cast<char*>(&rhs), tar, sizeof(rhs));
      galois_region_xor(reinterpret_cast<char*>(&lhs), tar, sizeof(rhs));
      m_sumCache[lhs][rhs] = m_sumCache[rhs][lhs];
    }
  }

  #ifdef DEBUG_GALOIS
    std::cerr<<"Galois: Generated GF sums table: "<<std::endl;
    std::cerr<<std::setfill(' ');
    std::cerr<<"      ";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(4)<<rhs;
    }
    std::cerr<<std::endl<<"      ";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(4)<<"____";
    }
    std::cerr<<std::endl;
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(5)<<rhs<<"|";
      for (unsigned int lhs = 0; lhs < getMax(); ++lhs) {
        std::cerr<<std::setw(4)<<m_sumCache[rhs][lhs];
      }
      std::cerr<<std::endl;
    }
  #endif
}

void Galois::precomputeProducts()
{
  m_productCache.resize(getMax(), std::vector<unsigned int>(getMax(), 0));
  for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
    for (unsigned int lhs = rhs; lhs < getMax(); ++lhs) {
      m_productCache[lhs][rhs] = m_productCache[rhs][lhs] =
        static_cast<unsigned int>(galois_single_multiply(
          static_cast<int>(rhs), static_cast<int>(lhs), static_cast<int>(m_w)
        ));
    }
  }

  #ifdef DEBUG_GALOIS
    std::cerr<<"Galois: Generated GF products table: "<<std::endl;
    std::cerr<<std::setfill(' ');
    std::cerr<<"      ";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(4)<<rhs;
    }
    std::cerr<<std::endl<<"      ";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(4)<<"____";
    }
    std::cerr<<std::endl;
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(5)<<rhs<<"|";
      for (unsigned int lhs = 0; lhs < getMax(); ++lhs) {
        std::cerr<<std::setw(4)<<m_productCache[rhs][lhs];
      }
      std::cerr<<std::endl;
    }
  #endif
}

void Galois::precomputeDivision()
{
  m_divideCache.resize(getMax(), std::vector<unsigned int>(getMax(), 0));
  for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
    for (unsigned int lhs = 1; lhs < getMax(); ++lhs) {
      m_divideCache[rhs][lhs] = static_cast<unsigned int>(
        galois_single_divide(
          static_cast<int>(rhs), static_cast<int>(lhs), static_cast<int>(m_w)
        )
      );
    }
  }

  #ifdef DEBUG_GALOIS
    std::cerr<<"Galois: Generated GF division table: "<<std::endl;
    std::cerr<<std::setfill(' ');
    std::cerr<<"      ";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(4)<<rhs;
    }
    std::cerr<<std::endl<<"      ";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(4)<<"____";
    }
    std::cerr<<std::endl;
    for (unsigned int rhs = 1; rhs < getMax(); ++rhs) {
      std::cerr<<std::setw(5)<<rhs<<"|";
      for (unsigned int lhs = 0; lhs < getMax(); ++lhs) {
        std::cerr<<std::setw(4)<<m_divideCache[rhs][lhs];
      }
      std::cerr<<std::endl;
    }
  #endif
}

void Galois::precomputeNumberOfOnes()
{
  m_numberOfOnesCache.resize(getMax(), 0);
  for (unsigned int rhs = 0; rhs < getMax(); ++rhs) {
    int *bitmatrix = jerasure_matrix_to_bitmatrix(1, 1, static_cast<int>(getW()), reinterpret_cast<int*>(&rhs));
    for (unsigned int el = 0; el < getW()*getW(); ++el) {
      m_numberOfOnesCache[rhs] += static_cast<unsigned int>(bitmatrix[el]);
    }
    delete[] bitmatrix;
  }

  #ifdef DEBUG_GALOIS
    std::cerr<<"Galois: Generated GF number of ones table: "<<std::endl;
    std::cerr<<std::setfill(' ');
    std::cerr<<"      ";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs)
      std::cerr<<std::setw(4)<<rhs;
    std::cerr<<std::endl<<"      ";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs)
      std::cerr<<std::setw(4)<<"____";
    std::cerr<<std::endl<<"  1's|";
    for (unsigned int rhs = 0; rhs < getMax(); ++rhs)
      std::cerr<<std::setw(4)<<m_numberOfOnesCache[rhs];
    std::cerr<<std::endl;
#endif
}

Galois::Galois(const unsigned int w):
  m_w(w)
{
  if (w > 12) throw std::runtime_error("We only support GF(2^w), where w <= 12");
  precomputeSums();
  precomputeProducts();
  precomputeDivision();
  precomputeNumberOfOnes();
  #ifdef DEBUG_GALOIS
    comparePerformance();
  #endif
}

#ifdef DEBUG_GALOIS
void Galois::comparePerformance() const
{
  using namespace std::chrono;
  {
    long long int sum1 = 0, sum2 = 0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for (unsigned int i = 0; i < 1<<25; ++i) {
      sum1 += sum(i % getMax(), (123*i) % getMax());
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    for(unsigned int i = 0; i < 1<<25; ++i) {
      sum2 += sum_uncached(i % getMax(), (123*i) % getMax());
    }
    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    if (sum1 == sum2) {
      std::cerr<<"Cached time for sum: "<<duration_cast<duration<double>>(t2 - t1).count()<<std::endl;
      std::cerr<<"Uncached time for sum: "<<duration_cast<duration<double>>(t3 - t2).count()<<std::endl;
    }
  }
  {
    long long int prod1 = 0, prod2 = 0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for (unsigned int i = 0; i < 1<<25; ++i) {
      prod1 += product(i % getMax(), (123*i) % getMax());
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    for(unsigned int i = 0; i < 1<<25; ++i) {
      prod2 += product_uncached(i % getMax(), (123*i) % getMax());
    }
    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    if (prod1 == prod2) {
      std::cerr<<"Cached time for product: "<<duration_cast<duration<double>>(t2 - t1).count()<<std::endl;
      std::cerr<<"Uncached time for product: "<<duration_cast<duration<double>>(t3 - t2).count()<<std::endl;
    }
  }
  {
    long long int div1 = 0, div2 = 0;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for (unsigned int i = 0; i < 1<<25; ++i) {
      div1 += divide(i % getMax(), (123*i) % getMax()==0?1:(123*i) % getMax());
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    for(unsigned int i = 0; i < 1<<25; ++i) {
      div2 += divide_uncached(i % getMax(), (123*i) % getMax()==0?1:(123*i) % getMax());
    }
    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    if (div1 == div2) {
      std::cerr<<"Cached time for division: "<<duration_cast<duration<double>>(t2 - t1).count()<<std::endl;
      std::cerr<<"Uncached time for division: "<<duration_cast<duration<double>>(t3 - t2).count()<<std::endl;
    }
  }
}
#endif
