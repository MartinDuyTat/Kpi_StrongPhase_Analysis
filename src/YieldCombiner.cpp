// Martin Duy Tat 14th July 2021

#include<string>
#include<vector>
#include<fstream>
#include<algorithm>
#include<numeric>
#include"YieldCombiner.h"
#include"BinVector.h"
#include"TMath.h"

YieldCombiner::YieldCombiner(int Bins, const std::vector<std::string> &Files): m_Yield(BinVector<double>(true, Bins)), m_YieldError(BinVector<double>(true, Bins)) {
  for(auto Filename : Files) {
    std::ifstream Infile(Filename);
    for(auto Bin : m_Yield.GetBinNumbers()) {
      double Yield, YieldStatError, YieldSystError, FlavourCorrectionError;
      Infile >> Yield >> YieldStatError >> YieldSystError >> FlavourCorrectionError;
      m_Yield[Bin] += Yield;
      m_YieldError[Bin] = TMath::Sqrt(TMath::Power(m_YieldError[Bin], 2) + TMath::Power(YieldStatError, 2) + TMath::Power(YieldSystError, 2) + TMath::Power(FlavourCorrectionError, 2));
    }
    Infile.close();
  }
  NormalizeYields();
}

void YieldCombiner::SaveYields(const std::string &Filename) const {
  std::ofstream Outfile(Filename);
  for(int Bin = 1; Bin <= m_Yield.Size(); Bin++) {
    Outfile << Bin << " " << m_Yield[Bin] << " " << m_YieldError[Bin] << " " << m_Yield[-Bin] << " " << m_YieldError[-Bin] << "\n";
  }
  Outfile.close();
}

void YieldCombiner::NormalizeYields() {
  double Sum = std::accumulate(m_Yield.begin(), m_Yield.end(), 0.0);
  std::transform(m_Yield.begin(), m_Yield.end(), m_Yield.begin(), [&Sum] (double a) { return a/Sum; });
  BinVector<double> YieldNormalizedError = m_YieldError;
  for(auto Bin : m_Yield.GetBinNumbers()) {
    m_YieldError[Bin] = CalculateNormalizationError(Bin, Sum, m_Yield, YieldNormalizedError);
  }
}

double YieldCombiner::CalculateNormalizationError(int Bin, double Sum, const BinVector<double> &NormalizedYield, const BinVector<double> &Error) const {
  double Error2 = 0.0;
  for(auto i : NormalizedYield.GetBinNumbers()) {
    if(i == Bin) {
      Error2 += TMath::Power((1 - NormalizedYield[i])*Error[i], 2);
    } else {
      Error2 += TMath::Power(NormalizedYield[i]*Error[i], 2);
    }
  }
  return TMath::Sqrt(Error2)/Sum;
}
