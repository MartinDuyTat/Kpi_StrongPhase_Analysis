// Martin Duy Tat

#include<fstream>
#include"AnalyseTruthYield.h"
#include"BinVector.h"
#include"TreeWrapper.h"

AnalyseTruthYield::AnalyseTruthYield(TreeWrapper *Tree): Analyse(Tree), m_GeneratorYields(true, m_BinningScheme.GetNumberBins()) {
}

BinVector<double> AnalyseTruthYield::CalculateTruthYield(const BinVector<double> &SignalMCYield, const std::string &Filename) {
  for(int i = 0; i < m_Tree->GetEntries(); i++) {
    m_Tree->GetEntry(i);
    int BinNumber = DetermineGeneratorBinNumber();
    if(BinNumber == 0) {
      continue;
    }
    m_GeneratorYields[i]++;
  }
  std::ofstream Outfile(Filename);
  BinVector<double> Efficiencies(true, m_BinningScheme.GetNumberBins());
  for(int i = 1; i <= m_BinningScheme.GetNumberBins; i++) {
    Efficiencies[i] = m_SignalMCYield[i]/m_GeneratorYields[i];
    Outfile << Efficiencies[i] << " ";
  }
  Outfile.close();
  return Efficiencies;
}
