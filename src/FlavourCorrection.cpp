// Martin Duy Tat 9th July 2021

#include<string>
#include<stdexcept>
#include<utility>
#include<fstream>
#include"FlavourCorrection.h"
#include"TMath.h"

FlavourCorrection::FlavourCorrection(int Bins, const std::string &Mode, const std::string &StrongPhaseFilename, const std::string &HadronicFilename): m_Bins(Bins), m_Mode(Mode), m_ci(BinVector<double>(false, m_Bins)), m_si(BinVector<double>(false, m_Bins)), m_Ki(BinVector<double>(true, m_Bins)) {
  if(m_Mode != "KSKK" && m_Mode != "KLKK") {
    throw std::invalid_argument("Unknown tag mode");
  }
  std::ifstream StrongPhaseFile(StrongPhaseFilename);
  for(int i = 1; i <= m_Bins; i++) {
    StrongPhaseFile >> m_ci[i] >> m_si[i] >> m_Ki[i] >> m_Ki[-i];
  }
  StrongPhaseFile.close();
  std::ifstream HadronicFile(HadronicFilename);
  std::vector<std::string> Parameters{"rD", "deltaD", "R"};
  for(auto iter = Parameters.begin(); iter != Parameters.end(); iter++) {
    double Value, Error;
    HadronicFile >> Value >> Error;
    m_Hadronic.insert({*iter, Value});
    m_HadronicErrors.insert({*iter, Error});
  }
  std::vector<std::string> Correlations{"rD_deltaD", "deltaD_R", "R_rD"};
  for(auto iter = Correlations.begin(); iter != Correlations.end(); iter++) {
    double Corr;
    HadronicFile >> Corr;
    m_HadronicCorrelations.insert({*iter, Corr});
  }
  HadronicFile.close();
}

std::pair<double, double> FlavourCorrection::CalculateCorrection(int Bin) const {
  int Sign = m_Mode == "KSKK" ? +1 : -1;
  int BinSign = Bin > 0 ? +1 : -1;
  int AbsBin = TMath::Abs(Bin);
  double KiKbari = TMath::Sqrt(m_Ki[Bin]*m_Ki[-Bin]);
  double Interference = TMath::Cos(m_Hadronic.at("deltaD"))*m_ci[AbsBin] - TMath::Sin(m_Hadronic.at("deltaD"))*BinSign*m_si[AbsBin];
  double Correction = m_Ki[Bin]/(m_Ki[Bin] + TMath::Power(m_Hadronic.at("rD"), 2)*m_Ki[-Bin] - 2*Sign*m_Hadronic.at("rD")*m_Hadronic.at("R")*KiKbari*Interference);
  double Error_rD = -(Correction*Correction/m_Ki[Bin])*2*(m_Hadronic.at("rD")*m_Ki[-Bin] - 2*Sign*m_Hadronic.at("R")*KiKbari*Interference);
  double Error_deltaD = -(Correction*Correction/m_Ki[Bin])*2*Sign*m_Hadronic.at("rD")*m_Hadronic.at("R")*KiKbari*(TMath::Sin(m_Hadronic.at("deltaD"))*m_ci[AbsBin] + TMath::Cos(m_Hadronic.at("deltaD"))*BinSign*m_si[AbsBin]);
  double Error_R = (Correction*Correction/m_Ki[Bin])*2*Sign*m_Hadronic.at("rD")*KiKbari*Interference;
  double Error = TMath::Sqrt(TMath::Power(Error_rD*m_HadronicErrors.at("rD"), 2)
                           + TMath::Power(Error_deltaD*m_HadronicErrors.at("deltaD"), 2)
                           + TMath::Power(Error_R*m_HadronicErrors.at("R"), 2))
                           + 2*Error_rD*Error_deltaD*m_HadronicErrors.at("rD")*m_HadronicErrors.at("deltaD")*m_HadronicCorrelations.at("rD_deltaD")
                           + 2*Error_deltaD*Error_R*m_HadronicErrors.at("deltaD")*m_HadronicErrors.at("R")*m_HadronicCorrelations.at("deltaD_R")
                           + 2*Error_R*Error_rD*m_HadronicErrors.at("R")*m_HadronicErrors.at("rD")*m_HadronicCorrelations.at("R_rD");
  return std::pair<double, double>{Correction, Error};
}

void FlavourCorrection::SaveCorrections(const std::string &Filename) const {
  std::ofstream File(Filename);
  for(auto Bin : m_Ki.GetBinNumbers()) {
    std::pair<double, double> Correction = CalculateCorrection(Bin);
    File << Bin << " " << Correction.first << " " << Correction.second << "\n";
  }
}
