// Martin Duy Tat 13th May 2021

#include<fstream>
#include<iostream>
#include"DoubleTagYields.h"

DoubleTagYields::DoubleTagYields(const std::string &Filename, int NBins): m_Yield(BinVector<double>(true, NBins)), m_YieldError(BinVector<double>(true, NBins)) {
  std::ifstream Infile(Filename);
  std::string line;
  while(std::getline(Infile, line)) {
    std::stringstream ss(line);
    int i;
    double Yield, YieldError;
    ss >> i >> Yield >> YieldError;
    m_Yield[i] = Yield;
    m_YieldError[i] = YieldError;
  }
  Infile.close();
}

double DoubleTagYields::GetYield(int Bin) const {
  return m_Yield[Bin];
}

double DoubleTagYields::GetYieldError(int Bin) const {
  return m_YieldError[Bin];
}

void DoubleTagYields::PrintYields() const {
  for(int i = 1; i <= m_Yield.Size(); i++) {
    std::cout << i << " " << GetYield(i) << " " << GetYieldError(i) << "\n";
  }
  for(int i = -1; i >= -m_Yield.Size(); i--) {
    std::cout << i << " " << GetYield(i) << " " << GetYieldError(i) << "\n";
  }
}
