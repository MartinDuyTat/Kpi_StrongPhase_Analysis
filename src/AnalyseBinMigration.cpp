// Martin Duy Tat 19th May 2021

#include<fstream>
#include<string>
#include"AnalyseBinMigration.h"

AnalyseBinMigration::AnalyseBinMigration(TreeWrapper *Tree): Analyse(Tree), m_BinYields(TMatrixT<double>(2*m_BinningScheme.GetNumberBins(), 2*m_BinningScheme.GetNumberBins())) {
}

void AnalyseBinMigration::CalculateBinMigrationYields(const std::string &Filename) {
  for(int i = 0; i < m_Tree->GetEntries(); i++) {
    m_Tree->GetEntry(i);
    int GeneratorBinNumber = DetermineGeneratorBinNumber();
    if(GeneratorBinNumber == 0) {
      continue;
    }
    int ReconstructedBinNumber = DetermineReconstructedBinNumber();
    if(ReconstructedBinNumber == 0) {
      ReconstructedBinNumber = DetermineMappedReconstructedBinNumber();
    }
    char Region = DetermineEnergyRegion();
    if(Region != 'S') {
      continue;
    }
    m_BinYields[ArrayIndex(GeneratorBinNumber)][ArrayIndex(ReconstructedBinNumber)]++;
  }
  SaveResults(Filename);
}
void AnalyseBinMigration::SaveResults(const std::string &Filename) const {
  std::ofstream Outfile(Filename);
  for(int i = 0; i < m_BinYields.GetNrows(); i++) {
    for(int j = 0; j < m_BinYields.GetNcols(); j++) {
      Outfile << m_BinYields[i][j] << " ";
    }
    Outfile << "\n";
  }
  Outfile.close();
}

TMatrixT<double> AnalyseBinMigration::GetBinMigrationMatrix() const {
  TMatrixT<double> BinMigrationMatrix = m_BinYields;
  for(int i = 0; i < m_BinYields.GetNrows(); i++) {
    double sum = 0.0;
    for(int j = 0; j < m_BinYields.GetNcols(); j++) {
      sum += m_BinYields[i][j];
    }
    for(int j = 0; j < m_BinYields.GetNcols(); j++) {
      BinMigrationMatrix[i][j] = m_BinYields[i][j]/sum;
    }
  }
  return BinMigrationMatrix.Invert();
}
