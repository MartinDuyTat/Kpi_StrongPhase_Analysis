// Martin Duy Tat 19th May 2021

#include<fstream>
#include<string>
#include<stdexcept>
#include"TreeWrapper.h"
#include"TEntryList.h"
#include"TCut.h"
#include"TRandom3.h"
#include"TMath.h"

TreeWrapper::TreeWrapper(const std::string &Filename, const std::string &TreeName, const std::string &CutFile, const std::string &DataType, double MomentumSmearing): m_Chain(TreeName.c_str()), m_MomentumSmearing(MomentumSmearing) {
  std::ifstream DataFile(Filename);
  std::string line;
  while(std::getline(DataFile, line)) {
    m_Chain.Add(line.c_str());
  }
  DataFile.close();
  TCut Cuts;
  if(CutFile != "None") {
    std::ifstream Infile(CutFile);
    while(std::getline(Infile, line)) {
      Cuts = Cuts && TCut(line.c_str());
    }
    Infile.close();
  }
  m_Chain.Draw(">> elist", Cuts, "entrylist");
  m_elist = (TEntryList*)gDirectory->Get("elist");
  m_Chain.SetEntryList(m_elist);
  if(DataType != "TruthTuple") {
    if(TreeName.find("KSKK") != std::string::npos) {
      m_SignalMode = "KSKK";
    } else if(TreeName.find("KLKK") != std::string::npos) {
      m_SignalMode = "KLKK";
    } else {
      throw std::runtime_error("Unknown signal mode");
    }
    if(TreeName.find("KeNu") != std::string::npos) {
      m_TagMode = "KeNu";
    } else if(TreeName.find("Kpipipi") != std::string::npos) {
      m_TagMode = "Kpipipi";
    } else if(TreeName.find("Kpipi0") != std::string::npos) {
      m_TagMode = "Kpipi0";
    } else if(TreeName.find("Kpi") != std::string::npos) {
      m_TagMode = "Kpi";
    }
  }
  SetBranchAddresses(DataType);
}

void TreeWrapper::SetBranchAddresses(const std::string &DataType) {
  bool RecData = DataType == "Data" || DataType == "SignalMC" || DataType == "TopoAna";
  bool GenData = DataType == "TruthTuple" || DataType == "SignalMC" || DataType == "TopoAna";
  if(RecData) {
    m_Chain.SetBranchAddress("SignalKalmanFitSuccess", &m_RecKinematics.KalmanFitSuccess);
    m_Chain.SetBranchAddress("TagKCharge", &m_RecKinematics.TagKCharge);
    if(m_SignalMode == "KSKK") {
      m_Chain.SetBranchAddress("SignalMBC", &m_RecKinematics.SignalMBC);
    } else if(m_SignalMode == "KLKK") {
      m_Chain.SetBranchAddress("SignalMMiss2", &m_RecKinematics.SignalMMiss2);
    }
    if(m_TagMode != "KeNu") {
      m_Chain.SetBranchAddress("TagMBC", &m_RecKinematics.TagMBC);
    } else {
      m_Chain.SetBranchAddress("TagUMiss", &m_RecKinematics.TagUMiss);
    }
    if(m_SignalMode == "KSKK") {
      m_Chain.SetBranchAddress("SignalKSpx", &m_RecKinematics.K0_P[0]);
      m_Chain.SetBranchAddress("SignalKSpy", &m_RecKinematics.K0_P[1]);
      m_Chain.SetBranchAddress("SignalKSpz", &m_RecKinematics.K0_P[2]);
      m_Chain.SetBranchAddress("SignalKSenergy", &m_RecKinematics.K0_P[3]);
      m_Chain.SetBranchAddress("SignalKSpxKalmanFit", &m_RecKinematics.K0Kalman_P[0]);
      m_Chain.SetBranchAddress("SignalKSpyKalmanFit", &m_RecKinematics.K0Kalman_P[1]);
      m_Chain.SetBranchAddress("SignalKSpzKalmanFit", &m_RecKinematics.K0Kalman_P[2]);
      m_Chain.SetBranchAddress("SignalKSenergyKalmanFit", &m_RecKinematics.K0Kalman_P[3]);
    } else if (m_SignalMode == "KLKK") {
      m_Chain.SetBranchAddress("SignalKLpx", &m_RecKinematics.K0_P[0]);
      m_Chain.SetBranchAddress("SignalKLpy", &m_RecKinematics.K0_P[1]);
      m_Chain.SetBranchAddress("SignalKLpz", &m_RecKinematics.K0_P[2]);
      m_Chain.SetBranchAddress("SignalKLenergy", &m_RecKinematics.K0_P[3]);
      m_Chain.SetBranchAddress("SignalKLpxKalmanFit", &m_RecKinematics.K0Kalman_P[0]);
      m_Chain.SetBranchAddress("SignalKLpyKalmanFit", &m_RecKinematics.K0Kalman_P[1]);
      m_Chain.SetBranchAddress("SignalKLpzKalmanFit", &m_RecKinematics.K0Kalman_P[2]);
      m_Chain.SetBranchAddress("SignalKLenergyKalmanFit", &m_RecKinematics.K0Kalman_P[3]);
    }
    m_Chain.SetBranchAddress("SignalKPluspx", &m_RecKinematics.KPlus_P[0]);
    m_Chain.SetBranchAddress("SignalKPluspy", &m_RecKinematics.KPlus_P[1]);
    m_Chain.SetBranchAddress("SignalKPluspz", &m_RecKinematics.KPlus_P[2]);
    m_Chain.SetBranchAddress("SignalKPlusenergy", &m_RecKinematics.KPlus_P[3]);
    m_Chain.SetBranchAddress("SignalKMinuspx", &m_RecKinematics.KMinus_P[0]);
    m_Chain.SetBranchAddress("SignalKMinuspy", &m_RecKinematics.KMinus_P[1]);
    m_Chain.SetBranchAddress("SignalKMinuspz", &m_RecKinematics.KMinus_P[2]);
    m_Chain.SetBranchAddress("SignalKMinusenergy", &m_RecKinematics.KMinus_P[3]);
    m_Chain.SetBranchAddress("SignalKPluspxKalmanFit", &m_RecKinematics.KPlusKalman_P[0]);
    m_Chain.SetBranchAddress("SignalKPluspyKalmanFit", &m_RecKinematics.KPlusKalman_P[1]);
    m_Chain.SetBranchAddress("SignalKPluspzKalmanFit", &m_RecKinematics.KPlusKalman_P[2]);
    m_Chain.SetBranchAddress("SignalKPlusenergyKalmanFit", &m_RecKinematics.KPlusKalman_P[3]);
    m_Chain.SetBranchAddress("SignalKMinuspxKalmanFit", &m_RecKinematics.KMinusKalman_P[0]);
    m_Chain.SetBranchAddress("SignalKMinuspyKalmanFit", &m_RecKinematics.KMinusKalman_P[1]);
    m_Chain.SetBranchAddress("SignalKMinuspzKalmanFit", &m_RecKinematics.KMinusKalman_P[2]);
    m_Chain.SetBranchAddress("SignalKMinusenergyKalmanFit", &m_RecKinematics.KMinusKalman_P[3]);
  }
  if(GenData) {
    m_Chain.SetBranchAddress("NumberOfParticles", &m_GenKinematics.NumberParticles);
    m_Chain.SetBranchAddress("ParticleIDs", m_GenKinematics.ParticleIDs.data());
    m_Chain.SetBranchAddress("MotherIndex", m_GenKinematics.MotherIndex.data());
    m_Chain.SetBranchAddress("True_Px", m_GenKinematics.TruePx.data());
    m_Chain.SetBranchAddress("True_Py", m_GenKinematics.TruePy.data());
    m_Chain.SetBranchAddress("True_Pz", m_GenKinematics.TruePz.data());
    m_Chain.SetBranchAddress("True_Energy", m_GenKinematics.TrueEnergy.data());
  }
  if(DataType == "TopoAna") {
    m_Chain.SetBranchAddress("iDcyTr", &m_GenKinematics.iDcyTr);
  }
}

const ReconstructedKinematics& TreeWrapper::GetReconstructedKinematics() const {
  return m_RecKinematics;
}

const GeneratorKinematics& TreeWrapper::GetGeneratorKinematics() const {
  return m_GenKinematics;
}

int TreeWrapper::GetEntries() const {
  return m_elist->GetN();
}

void TreeWrapper::GetEntry(int i) {
  m_Chain.GetEntry(m_Chain.GetEntryNumber(i));
  if(m_MomentumSmearing != 0.0) {
    SmearMomenta();
  }
}

std::string TreeWrapper::GetSignalMode() const {
  return m_SignalMode;
}

std::string TreeWrapper::GetTagMode() const {
  return m_TagMode;
}

void TreeWrapper::SmearMomenta() {
  for(int i = 0; i < 3; i++) {
    double RandomKplus = gRandom->Gaus(0.0, m_MomentumSmearing);
    double RandomKminus = gRandom->Gaus(0.0, m_MomentumSmearing);
    double RandomK0 = gRandom->Gaus(0.0, m_MomentumSmearing);
    m_RecKinematics.KPlus_P[i] += RandomKplus;
    m_RecKinematics.KPlusKalman_P[i] += RandomKplus;
    m_RecKinematics.KMinus_P[i] += RandomKminus;
    m_RecKinematics.KMinusKalman_P[i] += RandomKminus;
    m_RecKinematics.K0_P[i] += RandomK0;
    m_RecKinematics.K0Kalman_P[i] += RandomK0;
  }
  m_RecKinematics.KPlus_P[3] = TMath::Sqrt(TMath::Power(m_RecKinematics.KPlus_P[0], 2) + TMath::Power(m_RecKinematics.KPlus_P[1], 2) + TMath::Power(m_RecKinematics.KPlus_P[2], 2) + TMath::Power(0.493677, 2));
  m_RecKinematics.KPlusKalman_P[3] = TMath::Sqrt(TMath::Power(m_RecKinematics.KPlusKalman_P[0], 2) + TMath::Power(m_RecKinematics.KPlusKalman_P[1], 2) + TMath::Power(m_RecKinematics.KPlusKalman_P[2], 2) + TMath::Power(0.493677, 2));
  m_RecKinematics.KMinus_P[3] = TMath::Sqrt(TMath::Power(m_RecKinematics.KMinus_P[0], 2) + TMath::Power(m_RecKinematics.KMinus_P[1], 2) + TMath::Power(m_RecKinematics.KMinus_P[2], 2) + TMath::Power(0.493677, 2));
  m_RecKinematics.KMinusKalman_P[3] = TMath::Sqrt(TMath::Power(m_RecKinematics.KMinusKalman_P[0], 2) + TMath::Power(m_RecKinematics.KMinusKalman_P[1], 2) + TMath::Power(m_RecKinematics.KMinusKalman_P[2], 2) + TMath::Power(0.493677, 2));
  m_RecKinematics.K0_P[3] = TMath::Sqrt(TMath::Power(m_RecKinematics.K0_P[0], 2) + TMath::Power(m_RecKinematics.K0_P[1], 2) + TMath::Power(m_RecKinematics.K0_P[2], 2) + TMath::Power(0.493677, 2));
  m_RecKinematics.K0Kalman_P[3] = TMath::Sqrt(TMath::Power(m_RecKinematics.K0Kalman_P[0], 2) + TMath::Power(m_RecKinematics.K0Kalman_P[1], 2) + TMath::Power(m_RecKinematics.K0Kalman_P[2], 2) + TMath::Power(0.497611, 2));
}
