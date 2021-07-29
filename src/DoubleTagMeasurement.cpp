// Martin Duy Tat 13th May 2021

#include<iostream>
#include<string>
#include<vector>
#include"BinVector.h"
#include"DoubleTagMeasurement.h"

DoubleTagMeasurement::DoubleTagMeasurement(int NBins, const std::string &K0Mode, const std::string &cisiHadronicParametersFilename, const std::string &KiHadronicParametersFilename, const std::string &DTYieldsFilename): m_NBins(NBins), m_HParameters(cisiHadronicParametersFilename, KiHadronicParametersFilename, m_NBins, K0Mode), m_DTYields(DTYieldsFilename, m_NBins), m_Mode(K0Mode) {
  std::cout << "Adding the following strong phase parameters from the file " << cisiHadronicParametersFilename << "and fractional yields from the file " << KiHadronicParametersFilename << ":\n";
  m_HParameters.PrintHadronicParameters();
  std::cout << "Adding the following double tag yields from the file " << DTYieldsFilename << ":\n";
  m_DTYields.PrintYields();
}

double DoubleTagMeasurement::GetChi2(double Normalization, double rDcosDelta, double rDsinDelta, const std::string &ErrorCategory, const std::vector<std::pair<std::string, int>> &VetoBins) {
  BinVector<double> YieldPredictions, YieldKiErrorPredictions, YieldcisiErrorPredictions;
  m_HParameters.CalculateNormalizedYields(Normalization, rDcosDelta, rDsinDelta, YieldPredictions, YieldKiErrorPredictions, YieldcisiErrorPredictions);
  double Chi2 = 0.0;
  for(int Bin = -m_NBins; Bin <= m_NBins; Bin++ ) {
    if(Bin == 0) {
      continue;
    }
    if(VetoBins.size() != 0 && std::find(VetoBins.begin(), VetoBins.end(), std::pair<std::string, int>({m_Mode, Bin})) != VetoBins.end()) {
      continue;
    }
    double ErrorSquared = 0.0;
    if(ErrorCategory.find("Kpi") != std::string::npos) {
      ErrorSquared += TMath::Power(m_DTYields.GetYieldError(Bin), 2);
    }
    if(ErrorCategory.find("Ki") != std::string::npos) {
      ErrorSquared += TMath::Power(YieldKiErrorPredictions[Bin], 2);
    }
    if(ErrorCategory.find("cisi") != std::string::npos) {
      ErrorSquared += TMath::Power(YieldcisiErrorPredictions[Bin], 2);
    }
    Chi2 += TMath::Power(YieldPredictions[Bin] - m_DTYields.GetYield(Bin), 2)/ErrorSquared;
  }
  return Chi2;
}

int DoubleTagMeasurement::GetNBins() const {
  return m_NBins;
}

void DoubleTagMeasurement::SmearKi() {
  m_HParameters.SmearKi();
}

void DoubleTagMeasurement::Smearcisi(const std::vector<double> &Smearing) {
  m_HParameters.Smearcisi(Smearing);
}

void DoubleTagMeasurement::RemoveSmearing() {
  m_HParameters.RemoveSmearing();
}

std::string DoubleTagMeasurement::Mode() const {
  return m_Mode;
}
