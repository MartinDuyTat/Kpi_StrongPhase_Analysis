// Martin Duy Tat 13th May 2021

#include<iostream>
#include<string>
#include<vector>
#include"BinVector.h"
#include"DoubleTagMeasurement.h"

DoubleTagMeasurement::DoubleTagMeasurement(int NBins, const std::string &K0Mode, const std::string &HadronicParametersFilename, const std::string &DTYieldsFilename): m_NBins(NBins), m_HParameters(HadronicParametersFilename, m_NBins, K0Mode), m_DTYields(DTYieldsFilename, m_NBins) {
  std::cout << "Adding the following hadronic parameters from the file " << HadronicParametersFilename << ":\n";
  m_HParameters.PrintHadronicParameters();
  std::cout << "Adding the following double tag yields from the file " << DTYieldsFilename << ":\n";
  m_DTYields.PrintYields();
  std::string hh = (NBins == 8 ? "pipi" : "KK");
  m_Mode = K0Mode + hh;
}

double DoubleTagMeasurement::GetChi2(double Normalization, double rDcosDelta, double rDsinDelta, const std::string &ErrorCategory) {
  BinVector<double> YieldPredictions, YieldKiErrorPredictions, YieldcisiErrorPredictions;
  m_HParameters.CalculateNormalizedYields(Normalization, rDcosDelta, rDsinDelta, YieldPredictions, YieldKiErrorPredictions, YieldcisiErrorPredictions);
  double Chi2 = 0.0;
  for(int Bin = -m_NBins; Bin <= m_NBins; Bin++ ) {
    if(Bin == 0) {
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
