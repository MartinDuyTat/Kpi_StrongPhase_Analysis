// Martin Duy Tat 19th May 2021

#include<string>
#include<vector>
#include<utility>
#include<algorithm>
#include<fstream>
#include<sstream>
#include"BinVector.h"
#include"Settings.h"
#include"AnalyseYield.h"
#include"TreeWrapper.h"
#include"Analyse.h"
#include"TMatrix.h"
#include"TMath.h"
#include"TGraph.h"
#include"TCanvas.h"

AnalyseYield::AnalyseYield(TreeWrapper *Tree, bool SubtractBackground, const std::string &PeakingBackgroundFile): Analyse(Tree), m_DalitzCoordinates(BinVector<std::vector<std::pair<double, double>>>(true, m_BinningScheme.GetNumberBins())), m_SubtractBackground(SubtractBackground), m_EventsOutsideMBCSpace(0), m_EventsOutsidePhaseSpace(0) {
  m_Yields.insert({'S', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
  m_PeakingBackground.insert({'S', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
  if(m_Tree->GetSignalMode() == "KSKK" && m_Tree->GetTagMode() != "KeNu") {
    m_Yields.insert({'A', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
    m_Yields.insert({'B', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
    m_Yields.insert({'C', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
    m_Yields.insert({'D', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
  } else {
    m_Yields.insert({'L', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
    m_Yields.insert({'H', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
    m_PeakingBackground.insert({'L', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
    m_PeakingBackground.insert({'H', BinVector<double>(true, m_BinningScheme.GetNumberBins())});
  }
  if(PeakingBackgroundFile != "") {
    std::ifstream PeakingFile(PeakingBackgroundFile);
    std::string line;
    while(std::getline(PeakingFile, line)) {
      std::string iDcyTr;
      char Region;
      double Background;
      std::stringstream ss(line);
      ss >> iDcyTr >> Region;
      if((m_Tree->GetSignalMode() == "KSKK" && m_Tree->GetTagMode() != "KeNu") && iDcyTr == "Other") {
	continue;
      }
      for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
	ss >> Background;
	m_PeakingBackground.at(Region)[i] += Background;
      }
      for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
	ss >> Background;
	m_PeakingBackground.at(Region)[-i] += Background;
      }
    }
    PeakingFile.close();
  }
}

double AnalyseYield::GetBackgroundSubtractedYield(int Bin) const {
  if(m_Tree->GetSignalMode() == "KSKK" && m_Tree->GetTagMode() != "KeNu") {
    // Area in beam constrained 2D plane
    double A_S = 10*10;
    double A_A = 10*25;
    double A_B = 25*10;
    double A_D = 19.5*19.5;
    double A_C = 25*25 - 21.5*21.5;
    double Background = (A_S/A_D)*m_Yields.at('D')[Bin] + (A_S/A_A)*(m_Yields.at('A')[Bin] - (A_S/A_A)*m_Yields.at('D')[Bin]) + (A_S/A_B)*(m_Yields.at('B')[Bin] - (A_S/A_B)*m_Yields.at('D')[Bin]) + (A_S/A_C)*(m_Yields.at('C')[Bin] - (A_S/A_C)*m_Yields.at('D')[Bin]);
    return m_Yields.at('S')[Bin] - Background - m_PeakingBackground.at('S')[Bin];
  } else {
    double alpha = KpiSettings::Get().GetDouble("alpha");
    double beta = KpiSettings::Get().GetDouble("beta");
    double gamma = KpiSettings::Get().GetDouble("gamma");
    double delta = KpiSettings::Get().GetDouble("delta");
    return ((m_Yields.at('S')[Bin] - m_PeakingBackground.at('S')[Bin]) - delta*(m_Yields.at('L')[Bin] - m_PeakingBackground.at('L')[Bin]) - gamma*(m_Yields.at('H')[Bin] - m_PeakingBackground.at('H')[Bin]))/(1 - delta*alpha - gamma*beta);
  } 
}

void AnalyseYield::CalculateDoubleTagYields(const TMatrixT<double> &BinMigrationMatrix, const std::string &Filename) {
  for(int i = 0; i < m_Tree->GetEntries(); i++) {
    m_Tree->GetEntry(i);
    char Region = DetermineEnergyRegion();
    if(Region == 'F') {
      m_EventsOutsideMBCSpace++;
      continue;
    }
    int BinNumber = DetermineReconstructedBinNumber();
    if(BinNumber == 0) {
      m_EventsOutsidePhaseSpace++;
      BinNumber = DetermineMappedReconstructedBinNumber();
    }
    m_Yields.at(Region)[BinNumber]++;
    if(Region == 'S') {
      ReconstructedKinematics RecKinematics = m_Tree->GetReconstructedKinematics();
      double M2Plus, M2Minus;
      if(RecKinematics.KalmanFitSuccess == 1) {
	M2Plus = (RecKinematics.K0Kalman_P + RecKinematics.KPlusKalman_P).M2();
	M2Minus = (RecKinematics.K0Kalman_P + RecKinematics.KMinusKalman_P).M2();
      } else {
	M2Plus = (RecKinematics.K0_P + RecKinematics.KPlus_P).M2();
	M2Minus = (RecKinematics.K0_P + RecKinematics.KMinus_P).M2();
      }
      if(RecKinematics.TagKCharge < 0) {
	std::swap(M2Plus, M2Minus);
      }
      m_DalitzCoordinates[BinNumber].push_back(std::pair<double, double>{M2Plus, M2Minus});
    }
  }
  SaveFinalYields(BinMigrationMatrix, Filename);
}

void AnalyseYield::SaveFinalYields(const TMatrixT<double> &BinMigrationMatrix, const std::string &Filename) const {
  TMatrixT<double> RawDoubleTagYield(2*m_BinningScheme.GetNumberBins(), 1);
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    if(m_SubtractBackground) {
      RawDoubleTagYield[ArrayIndex(i)][0] = GetBackgroundSubtractedYield(i);
    } else {
      RawDoubleTagYield[ArrayIndex(i)][0] = m_Yields.at('S')[i];
    }
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    if(m_SubtractBackground) {
      RawDoubleTagYield[ArrayIndex(-i)][0] = GetBackgroundSubtractedYield(-i);
    } else {
      RawDoubleTagYield[ArrayIndex(-i)][0] = m_Yields.at('S')[-i];
    }
  }
  TMatrixT<double> MigrationCorrectedYield = BinMigrationMatrix*RawDoubleTagYield;
  std::ofstream ResultsFile(Filename);
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    ResultsFile << MigrationCorrectedYield[ArrayIndex(i)][0] << " ";
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    ResultsFile << MigrationCorrectedYield[ArrayIndex(-i)][0] << " ";
  }
  ResultsFile << "\n";
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    ResultsFile << RawDoubleTagYield[ArrayIndex(i)][0] << " ";
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    ResultsFile << RawDoubleTagYield[ArrayIndex(-i)][0] << " ";
  }
  ResultsFile << "\n" << m_EventsOutsideMBCSpace << " " << m_EventsOutsidePhaseSpace << "\n";
  std::vector<char> Regions;
  if(m_Tree->GetSignalMode() == "KSKK" && m_Tree->GetTagMode() != "KeNu") {
    Regions = std::vector<char>{'S', 'A', 'B', 'C', 'D'};
  } else {
    Regions = std::vector<char>{'S', 'L', 'H'};
  }
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

void AnalyseYield::SaveDalitzDistributions(const std::string &Filename) const {
  std::vector<int> BinList;
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    BinList.push_back(i);
  }
  for(int i = 1; i <= m_BinningScheme.GetNumberBins(); i++) {
    BinList.push_back(-i);
  }
  for(auto iter = BinList.begin(); iter != BinList.end(); iter++) {
    std::string DalitzFilename = Filename + std::string(*iter > 0 ? "p" : "m") + std::to_string(TMath::Abs(*iter)) + std::string(".png");
    std::vector<double> M2Plus, M2Minus;
    for(auto Dalitz_position : m_DalitzCoordinates[*iter]) {
      M2Plus.push_back(Dalitz_position.first);
      M2Minus.push_back(Dalitz_position.second);
    }
    TCanvas c("c", "c", 900, 900);
    m_BinningScheme.Draw("col");
    TGraph gr(M2Plus.size(), M2Plus.data(), M2Minus.data());
    gr.Draw("SAME *");
    c.SaveAs(DalitzFilename.c_str());
  }
}
