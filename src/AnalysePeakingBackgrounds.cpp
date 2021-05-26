// Martin Duy Tat 26th May 2021

#include<map>
#include<string>
#include<fstream>
#include"AnalysePeakingBackgrounds.h"
#include"BinVector.h"
#include"TreeWrapper.h"

AnalysePeakingBackgrounds::AnalysePeakingBackgrounds(TreeWrapper *Tree, const std::string &Filename): Analyse(Tree) {
  std::ifstream iDcyTrFile(Filename);
  int iDcyTr;
  while(iDcyTrFile >> iDcyTr) {
    m_PeakingBackgrounds.insert({iDcyTr, BinVector<int>(true, m_BinningScheme.GetNumberBins())});
  }
  iDcyTrFile.close();
}
void AnalysePeakingBackgrounds::CalculatePeakingBackgrounds(const std::string &Filename) {
  for(int i = 0; i < m_Tree->GetEntries(); i++) {
    m_Tree->GetEntry(i);
    int iDcyTr = m_Tree->GetGeneratorKinematics().iDcyTr;
    if(m_PeakingBackgrounds.find(iDcyTr) != m_PeakingBackgrounds.end()) {
      if(DetermineMBCRegion() != 'S') {
	continue;
      }
      int Bin = DetermineReconstructedBinNumber();
      if(Bin == 0) {
	continue;
      }
      m_PeakingBackgrounds.at(iDcyTr)[Bin]++;
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
  Outfile.close();
}
