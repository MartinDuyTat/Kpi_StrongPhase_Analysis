// Martin Duy Tat 19th May 2021

#include<fstream>
#include<algorithm>
#include"AnalyseFinalYields.h"
#include"BinVector.h"

AnalyseFinalYields::AnalyseFinalYields(const std::string &GeneratorYieldsFilename, const std::string &SignalMCYieldsFilename, const std::string &DataYieldsFilename): Analyse(nullptr), m_GeneratorYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_SignalMCYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_DataYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())) {
  std::ifstream GeneratorFile(GeneratorYieldsFilename);
  std::ifstream SignalMCFile(SignalMCYieldsFilename);
  std::ifstream DataFile(DataYieldsFilename);
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    GeneratorFile >> m_GeneratorYields[i];
    SignalMCFile >> m_SignalMCYields[i];
    DataFile >> m_DataYields[i];
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    GeneratorFile >> m_GeneratorYields[-i];
    SignalMCFile >> m_SignalMCYields[-i];
    DataFile >> m_DataYields[-i];
  }
}

void AnalyseFinalYields::CalculateFinalYields(const std::string &Filename) const {
  BinVector<double> FinalYields(true, m_BinningScheme.GetNumberBins());
  double Sum = 0.0;
  std::ofstream Outfile(Filename);
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    FinalYields[i] = m_DataYields[i]*m_GeneratorYields[i]/m_SignalMCYields[i];
    Sum += FinalYields[i];
    Outfile << FinalYields[i] << " ";
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    FinalYields[-i] = m_DataYields[-i]*m_GeneratorYields[-i]/m_SignalMCYields[-i];
    Sum += FinalYields[-i];
    Outfile << FinalYields[-i] << " ";
  }
  Outfile << "\n";
  std::transform(FinalYields.begin(), FinalYields.end(), FinalYields.begin(), [&Sum](double value){return value/Sum;});
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    Outfile << FinalYields[i] << " ";
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    Outfile << FinalYields[-i] << " ";
  }
  Outfile.close();
}
