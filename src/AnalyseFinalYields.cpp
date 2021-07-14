// Martin Duy Tat 19th May 2021

#include<fstream>
#include<sstream>
#include<algorithm>
#include"AnalyseFinalYields.h"
#include"BinVector.h"

AnalyseFinalYields::AnalyseFinalYields(const std::string &GeneratorYieldsFilename, const std::string &SignalMCYieldsFilename, const std::string &DataYieldsFilename, const std::string &FlavourTagCorrectionFilename): Analyse(nullptr), m_GeneratorYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_SignalMCYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_DataYields(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_DataYieldStatErrors(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_DataYieldSystErrors(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_FlavourTagCorrections(BinVector<double>(true, m_BinningScheme.GetNumberBins())), m_FlavourTagCorrectionErrors(BinVector<double>(true, m_BinningScheme.GetNumberBins())) {
  std::ifstream GeneratorFile(GeneratorYieldsFilename);
  std::ifstream SignalMCFile(SignalMCYieldsFilename);
  std::ifstream DataFile(DataYieldsFilename);
  for(int i : m_DataYields.GetBinNumbers()) {
    double dummy;
    GeneratorFile >> m_GeneratorYields[i];
    SignalMCFile >> m_SignalMCYields[i] >> dummy >> dummy;
    DataFile >> m_DataYields[i] >> m_DataYieldStatErrors[i] >> m_DataYieldSystErrors[i];
  }
  if(FlavourTagCorrectionFilename != "") {
    std::ifstream FlavourCorrectionFile(FlavourTagCorrectionFilename);
    std::string line;
    while(std::getline(FlavourCorrectionFile, line)) {
      int Bin;
      double Correction, Error;
      std::stringstream ss(line);
      ss >> Bin >> Correction >> Error;
      m_FlavourTagCorrections[Bin] = Correction;
      m_FlavourTagCorrectionErrors[Bin] = Error;
    }
    FlavourCorrectionFile.close();
  } else {
    for(auto i : m_FlavourTagCorrections.GetBinNumbers()) {
      m_FlavourTagCorrections[i] = 1.0;
      m_FlavourTagCorrectionErrors[i] = 0.0;
    }
  }
}

void AnalyseFinalYields::CalculateFinalYields(const std::string &Filename) const {
  BinVector<double> FinalYields(true, m_BinningScheme.GetNumberBins());
  BinVector<double> FinalYieldStatErrors(true, m_BinningScheme.GetNumberBins());
  BinVector<double> FinalYieldSystErrors(true, m_BinningScheme.GetNumberBins());
  double Sum = 0.0;
  std::ofstream Outfile(Filename);
  for(int i : FinalYields.GetBinNumbers()) {
    double p = m_SignalMCYields[i]/m_GeneratorYields[i];
    double pError = TMath::Sqrt(p*(1 - p)/m_GeneratorYields[i]);
    FinalYields[i] = m_DataYields[i]*m_FlavourTagCorrections[i]/p;
    FinalYieldStatErrors[i] = (m_DataYieldStatErrors[i]/m_DataYields[i])*FinalYields[i];
    FinalYieldSystErrors[i] = TMath::Sqrt(TMath::Power(m_DataYieldSystErrors[i]/m_DataYields[i], 2) + TMath::Power(pError/p, 2) + TMath::Power(m_FlavourTagCorrectionErrors[i]/m_FlavourTagCorrections[i], 2))*FinalYields[i];
    Sum += FinalYields[i];
    Outfile << FinalYields[i] << " " << FinalYieldStatErrors[i] << " " << FinalYieldSystErrors[i] << " ";
  }
  Outfile << "\n";
  std::transform(FinalYields.begin(), FinalYields.end(), FinalYields.begin(), [&Sum](double value){return value/Sum;});
  for(int i : FinalYields.GetBinNumbers()) {
    Outfile << FinalYields[i] << " " << CalculateNormalizationError(i, Sum, FinalYields, FinalYieldStatErrors) << " " << CalculateNormalizationError(i, Sum, FinalYields, FinalYieldSystErrors) << " ";
  }
  Outfile << "\n";
  Outfile.close();
}

double AnalyseFinalYields::CalculateNormalizationError(int Bin, double Sum, const BinVector<double> &NormalizedYield, const BinVector<double> &Error) const {
  double Error2 = 0.0;
  for(auto i : NormalizedYield.GetBinNumbers()) {
    if(i == Bin) {
      Error2 += TMath::Power((1 - NormalizedYield[i])*Error[i], 2);
    } else {
      Error2 += TMath::Power(NormalizedYield[i]*Error[i], 2);
    }
  }
  return TMath::Sqrt(Error2)/Sum;
}
