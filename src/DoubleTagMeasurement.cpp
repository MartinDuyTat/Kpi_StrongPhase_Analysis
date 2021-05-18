// Martin Duy Tat 13th May 2021

#include<iostream>
#include<string>
#include"BinVector.h"
#include"DoubleTagMeasurement.h"

DoubleTagMeasurement::DoubleTagMeasurement(int NBins, const std::string &K0Mode, const std::string &HadronicParametersFilename, const std::string &DTYieldsFilename): m_NBins(NBins), m_HParameters(HadronicParametersFilename, m_NBins, K0Mode), m_DTYields(DTYieldsFilename, m_NBins) {
  std::cout << "Adding the following hadronic parameters from the file " << HadronicParametersFilename << ":\n";
  m_HParameters.PrintHadronicParameters();
  std::cout << "Adding the following double tag yields from the file " << DTYieldsFilename << ":\n";
  m_DTYields.PrintYields();
}

double DoubleTagMeasurement::GetChi2(double rDcosDelta, double rDsinDelta) {
  BinVector<double> YieldPredictions, YieldErrorPredictions;
  m_HParameters.CalculateNormalizedYields(rDcosDelta, rDsinDelta, YieldPredictions, YieldErrorPredictions);
  double Chi2 = 0.0;
  for(int Bin = -m_NBins; Bin <= m_NBins; Bin++ ) {
    if(Bin == 0) {
      continue;
    }
    double ErrorSquared = TMath::Power(m_DTYields.GetYieldError(Bin), 2) + TMath::Power(YieldErrorPredictions[Bin], 2);
    Chi2 += TMath::Power(YieldPredictions[Bin] - m_DTYields.GetYield(Bin), 2)/ErrorSquared;
  }
  return Chi2;
}

int DoubleTagMeasurement::GetNBins() const {
  return m_NBins;
}
