//
//  main.cpp
//  OptiECRS
//
//  Created by Mykyta Makovenko on 8/30/20.
//  Copyright © 2020 Texas A&M University. All rights reserved.
//

#include <iostream>
#include "galois/galois.h"
#include "matrix/extendedcauchymatrix.h"
#include "procedures/randomheuristic.h"
#include "procedures/optimizebranchandbound.h"
#include "procedures/optimizemixedinteger.h"
#include "procedures/optimizemixedintegeralt.h"
#include "gurobi_c++.h"

Galois getGFFromUser() {
  unsigned int w;
  std::cout << "\tGalois field w: ";
  std::cin >> w;
  return Galois(w);
}

ExtendedCauchyMatrix getInitialFromUser() {
  Galois GF = getGFFromUser();
  unsigned int rows, cols;
  std::cout << "\tNumber of rows (informational packets): ";
  std::cin >> rows;
  std::cout << "\tNumber of cols (parity packets): ";
  std::cin >> cols;
  return ExtendedCauchyMatrix(rows, cols, GF);
}

ExtendedCauchyMatrix runAlgorithm() {
  auto initial = getInitialFromUser();
  std::cout << "Initial:" << std::endl << initial << std::endl;
  double timelimit;
  std::cout << "\tTimelimit (seconds): ";
  std::cin >> timelimit;
  std::string algo;
  std::cout << "\tAlgorithm (BB, RH, IP or IPAlt): ";
  std::cin >> algo;
  if (algo == "BB") {
    auto algo = OptimizeBranchAndBound(initial);
    auto res = algo.run(timelimit);
    std::cout << "Time taken by algorithm: " << algo.getTimeInMS() << " seconds."<< std::endl;
    return res;
  }
  if (algo == "RH") {
    unsigned int N;
    std::cout << "\tIterations limit: ";
    std::cin >> N;
    auto algo = RandomHeuristic(initial, N);
    return algo.run(timelimit);
  }
  if (algo == "IP") {
    auto algo = OptimizeMixedInteger(initial);
    return algo.run(timelimit);
  }
  if (algo == "IPAlt") {
    auto algo = OptimizeMixedIntegerAlt(initial);
    return algo.run(timelimit);
  }
  else {
    std::cout << "\tUnknown algo." << std::endl;
    return initial;
  }
}

int main(int argc, const char * argv[]) {
  std::cout << runAlgorithm() << std::endl;
  return 0;
}
