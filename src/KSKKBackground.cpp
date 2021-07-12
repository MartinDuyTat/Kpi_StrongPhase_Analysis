// Martin Duy Tat 12th July 2021

#include<fstream>
#include<algorithm>
#include<numeric>
#include"KSKKBackground.h"
#include"Settings.h"
#include"TMath.h"

KSKKBackground::KSKKBackground(): m_Yields(true, 2), m_YieldErrors(true, 2) {
  if(KpiSettings::Get().GetString("Mode") != "KLKK") {
    return;
  }
  double rD, deltaD, rDError, deltaDError, R, RError;
  std::ifstream HadronicFile(KpiSettings::Get().GetString("HadronicParametersFile"));
  HadronicFile >> rD >> rDError >> deltaD >> deltaDError >> R >> RError;
  HadronicFile.close();
  BinVector<double> ci(false, 2), si(false, 2), ci_Error(false, 2), si_Error(false, 2), Ki(true, 2);
  std::ifstream StrongPhaseFile(KpiSettings::Get().GetString("KSKKStrongPhaseFile"));
  StrongPhaseFile >> ci[1] >> ci_Error[1] >> si[1] >> si_Error[1] >> Ki[1] >> Ki[-1];
  StrongPhaseFile >> ci[2] >> ci_Error[2] >> si[2] >> si_Error[2] >> Ki[2] >> Ki[-2];
  for(auto Bin : ci.GetBinNumbers()) {
    int AbsBin = TMath::Abs(Bin);
    int BinSign = Bin > 0 ? +1 : -1;
    double Interference = ci[AbsBin]*TMath::Cos(deltaD) - BinSign*si[AbsBin]*TMath::Sin(deltaD);
    double SqrtKiKi = TMath::Sqrt(Ki[Bin]*Ki[-Bin]);
    m_Yields[Bin] = Ki[Bin] + rD*rD*Ki[-Bin] - 2*rD*R*SqrtKiKi*Interference;
    double Error_rD = 2*rD*Ki[-Bin] - 2*R*SqrtKiKi*Interference;
    double Error_deltaD = 2*rD*R*(ci[AbsBin]*TMath::Sin(deltaD) + BinSign*si[AbsBin]*TMath::Cos(deltaD));
    double Error_R = -2*rD*SqrtKiKi*Interference;
    double Error_ci = -2*rD*R*SqrtKiKi*TMath::Cos(deltaD);
    double Error_si = 2*rD*R*SqrtKiKi*TMath::Sin(deltaD);
    m_YieldErrors[Bin] = TMath::Sqrt(TMath::Power(Error_rD*rDError, 2)
			      + TMath::Power(Error_deltaD*deltaDError, 2)
			      + TMath::Power(Error_R*RError, 2)
			      + TMath::Power(Error_ci*ci_Error[AbsBin], 2)
			      + TMath::Power(Error_si*si_Error[AbsBin], 2));
  }
  double TotalYield = KpiSettings::Get().GetDouble("KSKKYield");
  double Sum = std::accumulate(m_Yields.begin(), m_Yields.end(), 0.0);
  std::transform(m_Yields.begin(), m_Yields.end(), m_Yields.begin(), [&](auto& element) {return element*TotalYield/Sum;});
  std::transform(m_YieldErrors.begin(), m_YieldErrors.end(), m_YieldErrors.begin(), [&](auto& element) {return element*TotalYield/Sum;});
}




double KSKKBackground::GetBinYield(int Bin) const {
  return m_Yields[Bin];
}

double KSKKBackground::GetBinYieldError(int Bin) const {
  return m_YieldErrors[Bin];
}
