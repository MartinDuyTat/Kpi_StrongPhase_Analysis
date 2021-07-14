// Martin Duy Tat 14th July 2021

#include<string>
#include<vector>
#include<fstream>
#include<algorithm>
#include<numeric>
#include"YieldCombiner.h"
#include"BinVector.h"
#include"TMath.h"

YieldCombiner::YieldCombiner(int Bins, const std::vector<std::string> &Files): m_Yield(BinVector<double>(true, Bins)), m_YieldStatError(BinVector<double>(true, Bins)), m_YieldSystError(BinVector<double>(true, Bins)) {
  for(auto Filename : Files) {
    std::ifstream Infile(Filename);
    for(auto Bin : m_Yield.GetBinNumbers()) {
      double Yield, YieldStatError, YieldSystError;
      Infile >> Yield >> YieldStatError >> YieldSystError;
      m_Yield[Bin] += Yield;
      m_YieldStatError[Bin] = TMath::Sqrt(TMath::Power(m_YieldStatError[Bin], 2) + TMath::Power(YieldStatError, 2));
      m_YieldSystError[Bin] = TMath::Sqrt(TMath::Power(m_YieldSystError[Bin], 2) + TMath::Power(YieldSystError, 2));
    }
    Infile.close();
  }
  NormalizeYields();
}

void YieldCombiner::SaveYields(const std::string &Filename) const {
  std::ofstream Outfile(Filename);
  for(auto Bin : m_Yield.GetBinNumbers()) {
    Outfile << Bin << " " << m_Yield[Bin] << " " << m_YieldStatError[Bin] << " " << m_YieldSystError[Bin] << "\n";
  }
  Outfile.close();
}

void YieldCombiner::NormalizeYields() {
  double Sum = std::accumulate(m_Yield.begin(), m_Yield.end(), 0.0);
  std::transform(m_Yield.begin(), m_Yield.end(), m_Yield.begin(), [&Sum] (double a) { return a/Sum; });
  BinVector<double> YieldNormalizedStatError = m_YieldStatError;
  BinVector<double> YieldNormalizedSystError = m_YieldSystError;
  for(auto Bin : m_Yield.GetBinNumbers()) {
    m_YieldStatError[Bin] = CalculateNormalizationError(Bin, Sum, m_Yield, YieldNormalizedStatError);
    m_YieldSystError[Bin] = CalculateNormalizationError(Bin, Sum, m_Yield, YieldNormalizedSystError);
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
