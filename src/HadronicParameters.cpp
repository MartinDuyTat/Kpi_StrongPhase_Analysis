// Martin Duy Tat 13th May 2021

#include<iostream>
#include<stdexcept>
#include<algorithm>
#include"HadronicParameters.h"
#include"BinVector.h"
#include"TMath.h"

HadronicParameters::HadronicParameters(): m_NBins(0) {
}

HadronicParameters::HadronicParameters(const std::string &Filename, int NBins, const std::string &K0Mode): m_NBins(NBins), m_K0Mode(K0Mode), m_Ki(BinVector<double>(true, m_NBins)), m_KiError(BinVector<double>(true, m_NBins)), m_ci(BinVector<double>(false, m_NBins)), m_si(BinVector<double>(false, m_NBins)) {
  std::ifstream Infile(Filename);
  std::string line;
  while(std::getline(Infile, line)) {
    std::stringstream ss(line);
    int i;
    double ci, si, Ki, KiError, Kbari, KbariError;
    ss >> i >> ci >> si >> Ki >> KiError >> Kbari >> KbariError;
    m_ci[i] = ci;
    m_si[i] = si;
    m_Ki[i] = Ki;
    m_Ki[-i] = Kbari;
    m_KiError[i] = KiError;
    m_KiError[-i] = KbariError;
  }
  Infile.close();
}

std::string HadronicParameters::CreateYieldFormula(int Bin) const {
  std::string Ki = std::to_string(m_Ki[Bin]);
  std::string Kbari = std::to_string(m_Ki[-Bin]);
  std::string ci = std::to_string(m_ci[Bin > 0 ? Bin : -Bin]);
  std::string si = std::to_string(m_si[Bin > 0 ? Bin : -Bin]);
  std::string Formula = std::string("@0*(") + Ki + std::string(" + @1*@1*") + Kbari + std::string(" - 2*@1*sqrt(") + Ki + std::string("*") + Kbari + std::string(")*(") + ci + std::string("*@2 - ") + si + std::string("*@3))");
  return Formula;
}

double HadronicParameters::Getci(int Bin) const {
  return m_ci[TMath::Abs(Bin)];
}

double HadronicParameters::Getsi(int Bin) const {
  return m_si[TMath::Abs(Bin)];
}

double HadronicParameters::GetKi(int Bin) const {
  return m_Ki[Bin];
}

double HadronicParameters::GetKiError(int Bin) const {
  return m_KiError[Bin];
}

double HadronicParameters::NormalizeYield(double rDcosDelta, double rDsinDelta) const {
  double YieldSum = 0.0;
  for(int i = 1; i <= m_NBins; i++) {
    YieldSum += CalculateYield(i, rDcosDelta, rDsinDelta);
    YieldSum += CalculateYield(-i, rDcosDelta, rDsinDelta);
  }
  return 1.0/YieldSum;
}

double HadronicParameters::CalculateYield(int Bin, double rDcosDelta, double rDsinDelta) const {
  if(m_K0Mode != "KS" && m_K0Mode != "KL") {
    throw std::invalid_argument("Invalid K0 mode");
  }
  int CP = m_K0Mode == "KS" ? +1 : -1;
  return m_Ki[Bin] + (rDcosDelta*rDcosDelta + rDsinDelta*rDsinDelta)*m_Ki[-Bin] - 2*CP*TMath::Sqrt(m_Ki[Bin]*m_Ki[-Bin])*(m_ci[TMath::Abs(Bin)]*rDcosDelta - m_si[TMath::Abs(Bin)]*rDsinDelta);
}

double HadronicParameters::CalculateYieldError(int Bin, double rDcosDelta, double rDsinDelta) const {
  double KiDerivative = 1 + TMath::Sqrt(m_Ki[-Bin]/m_Ki[Bin])*(m_ci[TMath::Abs(Bin)]*rDcosDelta - m_si[TMath::Abs(Bin)]*rDsinDelta);
  double KbariDerivative = (rDcosDelta*rDcosDelta + rDsinDelta*rDsinDelta) + TMath::Sqrt(m_Ki[Bin]/m_Ki[-Bin])*(m_ci[TMath::Abs(Bin)]*rDcosDelta - m_si[TMath::Abs(Bin)]*rDsinDelta);
  return TMath::Sqrt(TMath::Power(KiDerivative*m_KiError[Bin], 2) + TMath::Power(KbariDerivative*m_KiError[-Bin], 2));
}

void HadronicParameters::CalculateNormalizedYields(double Normalization, double rDcosDelta, double rDsinDelta, BinVector<double> &Yield, BinVector<double> &YieldError) const {
  Yield = BinVector<double>(true, m_NBins);
  YieldError = BinVector<double>(true, m_NBins);
  for(int i = 1; i <= m_NBins; i++) {
    Yield[i] = CalculateYield(i, rDcosDelta, rDsinDelta);
    Yield[-i] = CalculateYield(-i, rDcosDelta, rDsinDelta);
    YieldError[i] = CalculateYieldError(i, rDcosDelta, rDsinDelta);
    YieldError[-i] = CalculateYieldError(-i, rDcosDelta, rDsinDelta);
  }
  double NormalizationConstant = NormalizeYield(rDcosDelta, rDsinDelta);
  std::transform(Yield.begin(), Yield.end(), Yield.begin(), [&NormalizationConstant, Normalization](auto &c){return Normalization*NormalizationConstant*c;});
  std::transform(YieldError.begin(), YieldError.end(), YieldError.begin(), [&NormalizationConstant, &Normalization](auto &c){return Normalization*NormalizationConstant*c;});
}

void HadronicParameters::PrintHadronicParameters() const {
  for(int i = 1; i <= m_ci.Size(); i++) {
    std::cout << i << " " << Getci(i) << " " << Getsi(i) << " " << GetKi(i) << " " << GetKiError(i) << " " << GetKi(-i) << " " << GetKiError(-i) << "\n";
  }
}
