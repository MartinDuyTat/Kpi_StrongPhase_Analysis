// Martin Duy Tat 22nd July 2021

#include<string>
#include<vector>
#include<map>
#include<fstream>
#include<sstream>
#include<algorithm>
#include<stdexcept>
#include"cisiCovariance.h"
#include"TMatrixT.h"
#include"TDecompChol.h"
#include"TRandom3.h"

cisiCovariance::cisiCovariance() {
  m_Modes = std::map<std::string, bool>({{"KSKK", false}, {"KLKK", false}, {"KSpipi", false}, {"KLpipi", false}});
  m_StatCorr.insert({"K0KK", TMatrixT<double>(8, 8)});
  m_StatCorr.insert({"K0pipi", TMatrixT<double>(32, 32)});
  m_SystCorr.insert({"K0KK", TMatrixT<double>(8, 8)});
  m_SystCorr.insert({"K0pipi", TMatrixT<double>(32, 32)});
  m_StatUncertainties.insert({"K0KK", std::vector<double>(8)});
  m_StatUncertainties.insert({"K0pipi", std::vector<double>(32)});
  m_SystUncertainties.insert({"K0KK", std::vector<double>(8)});
  m_SystUncertainties.insert({"K0pipi", std::vector<double>(32)});
  m_Smearing.insert({"KSKK", std::vector<double>(4)});
  m_Smearing.insert({"KLKK", std::vector<double>(4)});
  m_Smearing.insert({"KSpipi", std::vector<double>(16)});
  m_Smearing.insert({"KLpipi", std::vector<double>(16)});
}

void cisiCovariance::AddDataset(const std::string &Mode, const std::string &HadronicParameterFilename) {
  m_Modes[Mode] = true;
  int NBins, Shift;
  if(Mode == "KSKK" || Mode == "KLKK") {
    NBins = 2;
    Shift = Mode == "KSKK" ? 0 : 4;
  } else {
    NBins = 8;
    Shift = Mode == "KSpipi" ? 0 : 16;
  }
  std::string K0hhMode = NBins == 2 ? "K0KK" : "K0pipi";
  std::ifstream Infile(HadronicParameterFilename);
  for(int i = 0; i < NBins; i++) {
    std::string line;
    std::getline(Infile, line);
    std::stringstream ss(line);
    double dummy, ci_Stat, ci_Syst, si_Stat, si_Syst;
    ss >> dummy >> dummy >> ci_Stat >> ci_Syst >> dummy >> si_Stat >> si_Syst;
    m_StatUncertainties[K0hhMode][i + Shift] = ci_Stat;
    m_StatUncertainties[K0hhMode][i + NBins + Shift] = si_Stat;
    m_SystUncertainties[K0hhMode][i + Shift] = ci_Syst;
    m_SystUncertainties[K0hhMode][i + NBins + Shift] = si_Syst;
  }
  if(Mode == "KLKK" || Mode == "KLpipi") {
    Infile.close();
    return;
  }
  for(int i = 0; i < 4*NBins; i++) {
    for(int j = 0; j < 4*NBins; j++) {
      if(i > j) {
	m_StatCorr[K0hhMode][i][j] = m_StatCorr[K0hhMode][j][i];
      } else {
	double Value;
	Infile >> Value;
	m_StatCorr[K0hhMode][i][j] = Value;
      }
    }
  }
  for(int i = 0; i < 4*NBins; i++) {
    for(int j = 0; j < 4*NBins; j++) {
      if(i > j) {
	m_SystCorr[K0hhMode][i][j] = m_SystCorr[K0hhMode][j][i];
      } else {
	double Value;
	Infile >> Value;
	m_SystCorr[K0hhMode][i][j] = Value;
      }
    }
  }
  Infile.close();
}

void cisiCovariance::PrepareCholesky() {
  if(!std::all_of(m_Modes.begin(), m_Modes.end(), [] (auto i) { return i.second; })) {
    throw std::runtime_error("Not all datasets loaded");
  }
  for(int i = 0; i < 32; i++) {
    for(int j = 0; j < 32; j++) {
      m_StatCorr["K0pipi"][i][j] = m_StatCorr["K0pipi"][i][j]*m_StatUncertainties["K0pipi"][i]*m_StatUncertainties["K0pipi"][j]/100.0;
      m_SystCorr["K0pipi"][i][j] = m_SystCorr["K0pipi"][i][j]*m_SystUncertainties["K0pipi"][i]*m_SystUncertainties["K0pipi"][j]/100.0;
    }
  }
  for(int i = 0; i < 8; i++) {
    for(int j = 0; j < 8; j++) {
      m_StatCorr["K0KK"][i][j] = m_StatCorr["K0KK"][i][j]*m_StatUncertainties["K0KK"][i]*m_StatUncertainties["K0KK"][j]/100.0;
      m_SystCorr["K0KK"][i][j] = m_SystCorr["K0KK"][i][j]*m_SystUncertainties["K0KK"][i]*m_SystUncertainties["K0KK"][j]/100.0;
    }
  }
  TDecompChol CholStatK0pipi(m_StatCorr["K0pipi"]);
  TDecompChol CholSystK0pipi(m_SystCorr["K0pipi"]);
  TDecompChol CholStatK0KK(m_StatCorr["K0KK"]);
  TDecompChol CholSystK0KK(m_SystCorr["K0KK"]);
  bool Success = CholStatK0pipi.Decompose();
  if(!Success) {
    throw std::runtime_error("K0pipi statistical covariance matrix not positive definite");
  }
  Success = CholSystK0pipi.Decompose();
  if(!Success) {
    throw std::runtime_error("K0pipi systematic covariance matrix not positive definite");
  }
  Success = CholStatK0KK.Decompose();
  if(!Success) {
    throw std::runtime_error("K0KK statistical covariance matrix not positive definite");
  }
  Success = CholSystK0KK.Decompose();
  if(!Success) {
    throw std::runtime_error("K0KK systematic covariance matrix not positive definite");
  }
  TMatrixT<double> Temp1 = CholStatK0pipi.GetU();
  Temp1 = Temp1.T();
  m_StatCholesky.insert({"K0pipi", Temp1});
  TMatrixT<double> Temp2 = CholStatK0KK.GetU();
  Temp2 = Temp2.T();
  m_StatCholesky.insert({"K0KK", Temp2});
  TMatrixT<double> Temp3 = CholSystK0pipi.GetU();
  Temp3 = Temp3.T();
  m_SystCholesky.insert({"K0pipi", Temp3});
  TMatrixT<double> Temp4 = CholSystK0KK.GetU();
  Temp4 = Temp4.T();
  m_SystCholesky.insert({"K0KK", Temp4});
}

void cisiCovariance::Smear() {
  TMatrixT<double> K0pipiStatSmearing(32, 1), K0pipiSystSmearing(32, 1), K0KKStatSmearing(8, 1), K0KKSystSmearing(8, 1);
  for(int i = 0; i < 32; i++) {
    K0pipiStatSmearing[i][0] = gRandom->Gaus();
    K0pipiSystSmearing[i][0] = gRandom->Gaus();
  }
  for(int i = 0; i < 8; i++) {
    K0KKStatSmearing[i][0] = gRandom->Gaus();
    K0KKSystSmearing[i][0] = gRandom->Gaus();
  }
  TMatrixT<double> K0pipiSmearing = m_StatCholesky["K0pipi"]*K0pipiStatSmearing + m_SystCholesky["K0pipi"]*K0pipiSystSmearing;
  TMatrixT<double> K0KKSmearing = m_StatCholesky["K0KK"]*K0KKStatSmearing + m_SystCholesky["K0KK"]*K0KKSystSmearing;
  for(int i = 0; i < 16; i++) {
    m_Smearing["KSpipi"][i] = K0pipiSmearing[i][0];
    m_Smearing["KLpipi"][i] = K0pipiSmearing[i + 16][0];
  }
  for(int i = 0; i < 4; i++) {
    m_Smearing["KSKK"][i] = K0KKSmearing[i][0];
    m_Smearing["KLKK"][i] = K0KKSmearing[i + 4][0];
  }
}

std::vector<double> cisiCovariance::GetSmearing(const std::string &Mode) const {
  return m_Smearing.at(Mode);
}
