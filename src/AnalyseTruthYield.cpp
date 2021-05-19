// Martin Duy Tat

#include<fstream>
#include"AnalyseTruthYield.h"
#include"BinVector.h"
#include"TreeWrapper.h"

AnalyseTruthYield::AnalyseTruthYield(TreeWrapper *Tree): Analyse(Tree), m_GeneratorYields(true, m_BinningScheme.GetNumberBins()) {
}

void AnalyseTruthYield::CalculateTruthYield(const std::string &Filename) {
  for(int i = 0; i < m_Tree->GetEntries(); i++) {
    m_Tree->GetEntry(i);
    int BinNumber = DetermineGeneratorBinNumber();
    if(BinNumber == 0) {
      continue;
    }
    m_GeneratorYields[BinNumber]++;
  }
  std::ofstream Outfile(Filename);
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    Outfile << m_GeneratorYields[i] << " ";
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    Outfile << m_GeneratorYields[-i] << " ";
  }
  Outfile.close();
}
