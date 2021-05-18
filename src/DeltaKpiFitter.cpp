// Martin Duy Tat 18th May 2021

#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include"DeltaKpiFitter.h"
#include"DoubleTagMeasurement.h"

DeltaKpiFitter::DeltaKpiFitter(const std::string &Filename) {
  std::ifstream Infile(Filename);
  std::string line;
  int N = 0;
  while(std::getline(Infile, line)) {
    if(line[0] == '#') {
      continue;
    }
    std::stringstream ss(line);
    int NBins;
    std::string K0Mode, HParameterFilename, DTYieldFilename;
    ss >> NBins >> K0Mode >> HParameterFilename >> DTYieldFilename;
    m_Measurements.AddMeasurement(DoubleTagMeasurement(NBins, K0Mode, HParameterFilename, DTYieldFilename));
    N++;
  }
  Infile.close();
  if(N == 0) {
    std::cout << "No datasets inputted\n";
  } else if(N == 1) {
    std::cout << "1 dataset added\n";
  } else {
    std::cout << N << " datasets added\n";
  }
}

void DeltaKpiFitter::RunFit(const std::string &Filename) {
  m_Measurements.MinimizeChi2();
  std::ofstream Outfile(Filename);
  Outfile << m_Measurements.GetChi2PerDegreesOfFreedom() << "\n";
  Outfile << m_Measurements.GetFittedrDcosDelta() << " " << m_Measurements.GetErrorrDcosDelta() << "\n";
  Outfile << m_Measurements.GetFittedrDsinDelta() << " " << m_Measurements.GetErrorrDsinDelta();
  Outfile.close();
}
