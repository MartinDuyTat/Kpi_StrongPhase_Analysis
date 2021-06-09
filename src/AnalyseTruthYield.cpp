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
    m_GeneratorYields[BinNumber]++;
  }
  std::ofstream Outfile(Filename);
  for(int i : m_GeneratorYields.GetBinNumbers()) {
    Outfile << m_GeneratorYields[i] << " ";
  }
  Outfile << "\n";
  Outfile.close();
}
