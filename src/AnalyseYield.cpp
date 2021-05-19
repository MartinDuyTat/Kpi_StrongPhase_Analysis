// Martin Duy Tat 19th May 2021

#include<string>
#include"AnalyseYield.h"
#include"Analyse.h"
#include"TMatrix.h"

AnalyseYield::AnalyseYield(TreeWrapper *Tree): Analyse(Tree), m_EventsOutsideMBCSpace(0), m_EventsOutsidePhaseSpace(0) {
  m_Yields.insert({'S', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
  m_Yields.insert({'A', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
  m_Yields.insert({'B', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
  m_Yields.insert({'C', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
  m_Yields.insert({'D', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
}

double AnalyseYield::GetBackgroundSubtractedYield(int Bin) const {
  // Area in beam constrained 2D plane
  double A_S = 10*10;
  double A_A = 10*25;
  double A_B = 25*10;
  double A_D = 19.5*19.5;
  double A_C = 25*25 - 21.5*21.5;
  double Background = (A_S/A_D)*m_Yields.at('D')[Bin] + (A_S/A_A)*(m_Yields.at('A')[Bin] - (A_S/A_A)*m_Yields.at('D')[Bin]) + (A_S/A_B)*(m_Yields.at('B')[Bin] - (A_S/A_B)*m_Yields.at('D')[Bin]) + (A_S/A_C)*(m_Yields.at('C')[Bin] - (A_S/A_C)*m_Yields.at('D')[Bin]);
  return m_Yields.at('S')[Bin] - Background;
}

void AnalyseYield::CalculateDoubleTagYields(const TMatrixT<double> &BinMigrationMatrix, const std::string &Filename) {
  for(int i = 0; i < m_Tree->GetEntries(); i++) {
    m_Tree->GetEntry(i);
    int BinNumber = DetermineReconstructedBinNumber();
    if(BinNumber == 0) {
      m_EventsOutsidePhaseSpace++;
      continue;
    }
    char Region = DetermineMBCRegion();
    if(Region == 'F') {
      m_EventsOutsideMBCSpace++;
      continue;
    }
    m_Yields.at(Region)[BinNumber]++;
  }
  SaveFinalYields(BinMigrationMatrix, Filename);
}

void AnalyseYield::SaveFinalYields(const TMatrixT<double> &BinMigrationMatrix, const std::string &Filename) const {
  TMatrixT<double> RawDoubleTagYield(2*m_BinningScheme.GetNumberBins(), 1);
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    RawDoubleTagYield[ArrayIndex(i)][0] = GetBackgroundSubtractedYield(i);
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    RawDoubleTagYield[ArrayIndex(-i)][0] = GetBackgroundSubtractedYield(-i);
  }
  RawDoubleTagYield = BinMigrationMatrix*RawDoubleTagYield;
  std::ofstream ResultsFile(Filename);
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    ResultsFile << RawDoubleTagYield[ArrayIndex(i)][0] << " ";
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    ResultsFile << RawDoubleTagYield[ArrayIndex(-i)][0] << " ";
  }
  ResultsFile << "\n";
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    ResultsFile << GetBackgroundSubtractedYield(i) << " ";
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    ResultsFile << GetBackgroundSubtractedYield(-i) << " ";
  }
  ResultsFile << "\n" << m_EventsOutsideMBCSpace << " " << m_EventsOutsidePhaseSpace << "\n";
  std::vector<char> Regions{'S', 'A', 'B', 'C', 'D'};
  for(auto R : Regions) {
    for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
      ResultsFile << m_Yields.at(R)[i] << " ";
    }
    for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
      ResultsFile << m_Yields.at(R)[-i] << " ";
    }
    ResultsFile << "\n";
  }
}
