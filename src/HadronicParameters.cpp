// Martin Duy Tat 13th May 2021

#include<iostream>
#include<stdexcept>
#include<algorithm>
#include"HadronicParameters.h"
#include"BinVector.h"
#include"TMath.h"
#include"TRandom3.h"

HadronicParameters::HadronicParameters(): m_NBins(0) {
}

HadronicParameters::HadronicParameters(const std::string &Filename, int NBins, const std::string &K0Mode): m_NBins(NBins), m_K0Mode(K0Mode.substr(0, 2)), m_Ki(BinVector<double>(true, m_NBins)), m_KiError(BinVector<double>(true, m_NBins)), m_ci(BinVector<double>(false, m_NBins)), m_ciError(BinVector<double>(false, m_NBins)), m_si(BinVector<double>(false, m_NBins)), m_siError(BinVector<double>(false, m_NBins)), m_KiSmearing(BinVector<double>(true, m_NBins)), m_ciSmearing(BinVector<double>(false, m_NBins)), m_siSmearing(BinVector<double>(false, m_NBins)), m_Covariance_Cholesky(TMatrixT<double>(2*m_NBins, 2*m_NBins)) {
  std::ifstream Infile(Filename);
  std::string line;
  for(int i = 1; i <= m_NBins; i++) {
    std::getline(Infile, line);
    std::stringstream ss(line);
    double ci, ciErrorStat, ciErrorSyst, si, siErrorStat, siErrorSyst, Ki, KiError, Kbari, KbariError, dummy;
    ss >> dummy >> ci >> ciErrorStat >> ciErrorSyst >> si >> siErrorStat >> siErrorSyst >> Ki >> KiError >> Kbari >> KbariError;
    m_ci[i] = ci;
    m_ciError[i] = TMath::Sqrt(ciErrorStat*ciErrorStat + ciErrorSyst*ciErrorSyst);
    m_si[i] = si;
    m_siError[i] = TMath::Sqrt(siErrorStat*siErrorStat + siErrorSyst*siErrorSyst);
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
  return m_ci[TMath::Abs(Bin)] + m_ciSmearing[TMath::Abs(Bin)];
}

double HadronicParameters::Getsi(int Bin) const {
  int Sign = Bin > 0 ? +1 : -1;
  return Sign*(m_si[TMath::Abs(Bin)] + m_siSmearing[TMath::Abs(Bin)]);
}

double HadronicParameters::GetKi(int Bin) const {
  return m_Ki[Bin] + m_KiSmearing[Bin];
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
  return GetKi(Bin) + (rDcosDelta*rDcosDelta + rDsinDelta*rDsinDelta)*GetKi(-Bin) - 2*CP*TMath::Sqrt(GetKi(Bin)*GetKi(-Bin))*(Getci(Bin)*rDcosDelta - Getsi(Bin)*rDsinDelta);
}

double HadronicParameters::CalculateYieldKiError(int Bin, double rDcosDelta, double rDsinDelta) const {
  int CP = m_K0Mode == "KS" ? +1 : -1;
  double KiDerivative = 1 - CP*TMath::Sqrt(GetKi(-Bin)/GetKi(Bin))*(Getci(Bin)*rDcosDelta - Getsi(Bin)*rDsinDelta);
  double KbariDerivative = (rDcosDelta*rDcosDelta + rDsinDelta*rDsinDelta) - CP*TMath::Sqrt(GetKi(Bin)/GetKi(-Bin))*(Getci(Bin)*rDcosDelta - Getsi(Bin)*rDsinDelta);
  return TMath::Sqrt(TMath::Power(KiDerivative*m_KiError[Bin], 2) + TMath::Power(KbariDerivative*m_KiError[-Bin], 2));
}

double HadronicParameters::CalculateYieldcisiError(int Bin, double rDcosDelta, double rDsinDelta) const {
  int CP = m_K0Mode == "KS" ? +1 : -1;
  double ciDerivative = -CP*2*TMath::Sqrt(GetKi(Bin)*GetKi(-Bin))*rDcosDelta;
  double siDerivative = CP*2*TMath::Sqrt(GetKi(Bin)*GetKi(-Bin))*rDsinDelta;
  return TMath::Sqrt(TMath::Power(ciDerivative*m_ciError[TMath::Abs(Bin)], 2) + TMath::Power(siDerivative*m_siError[TMath::Abs(Bin)], 2));
}

void HadronicParameters::CalculateNormalizedYields(double Normalization, double rDcosDelta, double rDsinDelta, BinVector<double> &Yield, BinVector<double> &YieldKiError, BinVector<double> &YieldcisiError) const {
  Yield = BinVector<double>(true, m_NBins);
  YieldKiError = BinVector<double>(true, m_NBins);
  YieldcisiError = BinVector<double>(true, m_NBins);
  for(auto Bin : Yield.GetBinNumbers()) {
    Yield[Bin] = CalculateYield(Bin, rDcosDelta, rDsinDelta);
    YieldKiError[Bin] = CalculateYieldKiError(Bin, rDcosDelta, rDsinDelta);
    YieldcisiError[Bin] = CalculateYieldcisiError(Bin, rDcosDelta, rDsinDelta);
  }
  double NormalizationConstant = NormalizeYield(rDcosDelta, rDsinDelta);
  std::transform(Yield.begin(), Yield.end(), Yield.begin(), [&NormalizationConstant, Normalization](auto &c){return Normalization*NormalizationConstant*c;});
  std::transform(YieldKiError.begin(), YieldKiError.end(), YieldKiError.begin(), [&NormalizationConstant, &Normalization](auto &c){return Normalization*NormalizationConstant*c;});
  std::transform(YieldcisiError.begin(), YieldcisiError.end(), YieldcisiError.begin(), [&NormalizationConstant, &Normalization](auto &c){return Normalization*NormalizationConstant*c;});
}

void HadronicParameters::PrintHadronicParameters() const {
  for(int i = 1; i <= m_ci.Size(); i++) {
    std::cout << i << " " << Getci(i) << " " << Getsi(i) << " " << GetKi(i) << " " << GetKiError(i) << " " << GetKi(-i) << " " << GetKiError(-i) << "\n";
  }
}

void HadronicParameters::SmearKi() {
  for(auto Bin : m_KiSmearing.GetBinNumbers()) {
    m_KiSmearing[Bin] = gRandom->Gaus(0.0, m_KiError[Bin]);
  }
}

void HadronicParameters::Smearcisi(const std::vector<double> &Smearing) {
  for(int i = 1; i <= m_NBins; i++) {
    m_ciSmearing[i] = Smearing[i - 1];
    m_siSmearing[i] = Smearing[i + m_NBins - 1];
  }
}

void HadronicParameters::RemoveSmearing() {
  std::fill(m_KiSmearing.begin(), m_KiSmearing.end(), 0.0);
  std::fill(m_ciSmearing.begin(), m_ciSmearing.end(), 0.0);
  std::fill(m_siSmearing.begin(), m_siSmearing.end(), 0.0);
}
