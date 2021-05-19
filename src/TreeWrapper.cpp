// Martin Duy Tat 19th May 2021

#include<fstream>
#include<string>
#include"TreeWrapper.h"
#include"TEntryList.h"
#include"TCut.h"

TreeWrapper::TreeWrapper(const std::string &Filename, const std::string &TreeName, const std::string &CutFile, const std::string &DataType): m_Chain(TreeName.c_str()) {
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
  SetBranchAddresses(DataType);
}

void TreeWrapper::SetBranchAddresses(const std::string &DataType) {
  bool RecData = DataType == "Data" || DataType == "SignalMC";
  bool GenData = DataType == "TruthTuple" || DataType == "SignalMC";
  if(RecData) {
    m_Chain.SetBranchAddress("SignalKalmanFitSuccess", &m_RecKinematics.KalmanFitSuccess);
    m_Chain.SetBranchAddress("TagKCharge", &m_RecKinematics.TagKCharge);
    m_Chain.SetBranchAddress("SignalMBC", &m_RecKinematics.SignalMBC);
    m_Chain.SetBranchAddress("TagMBC", &m_RecKinematics.TagMBC);
    m_Chain.SetBranchAddress("SignalKSpx", &m_RecKinematics.KS_P[0]);
    m_Chain.SetBranchAddress("SignalKSpy", &m_RecKinematics.KS_P[1]);
    m_Chain.SetBranchAddress("SignalKSpz", &m_RecKinematics.KS_P[2]);
    m_Chain.SetBranchAddress("SignalKSenergy", &m_RecKinematics.KS_P[3]);
    m_Chain.SetBranchAddress("SignalKPluspx", &m_RecKinematics.KPlus_P[0]);
    m_Chain.SetBranchAddress("SignalKPluspy", &m_RecKinematics.KPlus_P[1]);
    m_Chain.SetBranchAddress("SignalKPluspz", &m_RecKinematics.KPlus_P[2]);
    m_Chain.SetBranchAddress("SignalKPlusenergy", &m_RecKinematics.KPlus_P[3]);
    m_Chain.SetBranchAddress("SignalKMinuspx", &m_RecKinematics.KMinus_P[0]);
    m_Chain.SetBranchAddress("SignalKMinuspy", &m_RecKinematics.KMinus_P[1]);
    m_Chain.SetBranchAddress("SignalKMinuspz", &m_RecKinematics.KMinus_P[2]);
    m_Chain.SetBranchAddress("SignalKMinusenergy", &m_RecKinematics.KMinus_P[3]);
    m_Chain.SetBranchAddress("SignalKSpxKalmanFit", &m_RecKinematics.KSKalman_P[0]);
    m_Chain.SetBranchAddress("SignalKSpyKalmanFit", &m_RecKinematics.KSKalman_P[1]);
    m_Chain.SetBranchAddress("SignalKSpzKalmanFit", &m_RecKinematics.KSKalman_P[2]);
    m_Chain.SetBranchAddress("SignalKSenergyKalmanFit", &m_RecKinematics.KSKalman_P[3]);
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
}
