// Martin Duy Tat 26th May 2021

#include<map>
#include<string>
#include<fstream>
#include<sstream>
#include<algorithm>
#include"AnalysePeakingBackgrounds.h"
#include"BinVector.h"
#include"TreeWrapper.h"

AnalysePeakingBackgrounds::AnalysePeakingBackgrounds(TreeWrapper *Tree, const std::string &Filename): Analyse(Tree), m_OtherBackgrounds(BinVector<double>(true, m_BinningScheme.GetNumberBins())) {
  std::ifstream iDcyTrFile(Filename);
  std::string line1, line2, line3;
  std::getline(iDcyTrFile, line1);
  std::getline(iDcyTrFile, line2);
  std::getline(iDcyTrFile, line3);
  int iDcyTr;
  std::stringstream ss1(line1), ss2(line2), ss3(line3);
  while(ss1 >> iDcyTr) {
    m_SignalComponents.push_back(iDcyTr);
  }
  while(ss2 >> iDcyTr) {
    double ScaleFactor;
    ss3 >> ScaleFactor;
    m_PeakingBackgrounds.insert({iDcyTr, BinVector<double>(true, m_BinningScheme.GetNumberBins())});
    m_ScaleFactors.insert({iDcyTr, ScaleFactor});
  }
  iDcyTrFile.close();
}

void AnalysePeakingBackgrounds::CalculatePeakingBackgrounds(const std::string &Filename, double MCScale) {
  for(int i = 0; i < m_Tree->GetEntries(); i++) {
    m_Tree->GetEntry(i);
    int iDcyTr = m_Tree->GetGeneratorKinematics().iDcyTr;
    if(std::find(m_SignalComponents.begin(), m_SignalComponents.end(), iDcyTr) == m_SignalComponents.end()) {
      if(DetermineMBCRegion() != 'S') {
	continue;
      }
      int Bin = DetermineReconstructedBinNumber();
      if(Bin == 0) {
	continue;
      }
      if(m_PeakingBackgrounds.find(iDcyTr) != m_PeakingBackgrounds.end()) {
	m_PeakingBackgrounds.at(iDcyTr)[Bin] += m_ScaleFactors.at(iDcyTr)/MCScale;
      } else {
	m_OtherBackgrounds[Bin] += 1.0/MCScale;
      }
    }
  }
  SavePeakingBackgrounds(Filename);
}

void AnalysePeakingBackgrounds::SavePeakingBackgrounds(const std::string &Filename) const {
  std::ofstream Outfile(Filename);
  for(auto iter = m_PeakingBackgrounds.begin(); iter != m_PeakingBackgrounds.end(); iter++) {
    Outfile << iter->first << " ";
    for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
      Outfile << iter->second[i] << " ";
    }
    for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
      Outfile << iter->second[-i] << " ";
    }
    Outfile << "\n";
  }
  Outfile << "Other ";
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    Outfile << m_OtherBackgrounds[i] << " ";
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    Outfile << m_OtherBackgrounds[-i] << " ";
  }
  Outfile.close();
}
