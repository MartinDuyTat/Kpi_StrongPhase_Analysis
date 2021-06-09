// Martin Duy Tat 19th May 2021

#include<fstream>
#include<algorithm>
#include"AnalyseFinalYields.h"
#include"BinVector.h"

AnalyseFinalYields::AnalyseFinalYields(const std::string &GeneratorYieldsFilename, const std::string &SignalMCYieldsFilename, const std::string &DataYieldsFilename): Analyse(nullptr), m_GeneratorYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_SignalMCYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_DataYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_DataYieldErrors(BinVector<double>(true, m_BinningScheme.GetNumberBins())) {
  std::ifstream GeneratorFile(GeneratorYieldsFilename);
  std::ifstream SignalMCFile(SignalMCYieldsFilename);
  std::ifstream DataFile(DataYieldsFilename);
  for(int i : m_DataYields.GetBinNumbers()) {
    double dummy;
    GeneratorFile >> m_GeneratorYields[i];
    SignalMCFile >> m_SignalMCYields[i] >> dummy;
    DataFile >> m_DataYields[i] >> m_DataYieldErrors[i];
  }
}

void AnalyseFinalYields::CalculateFinalYields(const std::string &Filename) const {
  BinVector<double> FinalYields(true, m_BinningScheme.GetNumberBins());
  BinVector<double> FinalYieldErrors(true, m_BinningScheme.GetNumberBins());
  double Sum = 0.0;
  std::ofstream Outfile(Filename);
  for(int i : FinalYields.GetBinNumbers()) {
    double p = m_SignalMCYields[i]/m_GeneratorYields[i];
    double pError = TMath::Sqrt(p*(1 - p)/m_GeneratorYields[i]);
    FinalYields[i] = m_DataYields[i]/p;
    FinalYieldErrors[i] = TMath::Sqrt(TMath::Power(m_DataYieldErrors[i]/m_DataYields[i], 2) + TMath::Power(pError/p, 2))*FinalYields[i];
    Sum += FinalYields[i];
    Outfile << FinalYields[i] << " " << FinalYieldErrors[i] << " ";
  }
  Outfile << "\n";
  std::transform(FinalYields.begin(), FinalYields.end(), FinalYields.begin(), [&Sum](double value){return value/Sum;});
  for(int i : FinalYields.GetBinNumbers()) {
    Outfile << FinalYields[i] << " " << FinalYieldErrors[i]*Sum*(1 - FinalYields[i])/(Sum*Sum) << " ";
  }
  Outfile << "\n";
  Outfile.close();
}
