// Martin Duy Tat 19th May 2021
/**
 * AnalyseDoubleTagYields is an application for calculating the double tag yield and correct for bin migration and bin efficiencies
 * @param 1 Filename of settings file
 */

#include<iostream>
#include<string>
#include"Settings.h"
#include"TreeWrapper.h"
#include"AnalyseYield.h"
#include"AnalyseBinMigration.h"
#include"AnalyseTruthYield.h"
#include"AnalyseFinalYields.h"

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 1 input argument\n";
    return 0;
  }
  KpiSettings::Get().Initialize(std::string(argv[1]));
  TreeWrapper MCSignalTree(KpiSettings::Get().GetString("SignalMCFiles"), KpiSettings::Get().GetString("TreeName"), KpiSettings::Get().GetString("TruthMatchingCutFile"), "SignalMC", KpiSettings::Get().GetDouble("SignalMCMomentumSmearing"));
  AnalyseBinMigration AnalysisBinMigration(&MCSignalTree);
  AnalysisBinMigration.CalculateBinMigrationYields(KpiSettings::Get().GetString("BinMigrationMatrixFile"));
  AnalyseYield AnalysisSignalMC(&MCSignalTree, false, "", true);
  AnalysisSignalMC.CalculateDoubleTagYields(AnalysisBinMigration.GetBinMigrationMatrix(), KpiSettings::Get().GetString("SignalMCResultsFile"));
  TreeWrapper TruthTree(KpiSettings::Get().GetString("TruthTupleFiles"), KpiSettings::Get().GetString("TruthTreeName"), "None", "TruthTuple");
  AnalyseTruthYield AnalysisTruthTuple(&TruthTree);
  AnalysisTruthTuple.CalculateTruthYield(KpiSettings::Get().GetString("TruthResultsFile"));
  TreeWrapper DataTree(KpiSettings::Get().GetString("DataFiles"), KpiSettings::Get().GetString("TreeName"), KpiSettings::Get().GetString("DataCutFile"), "Data");
  AnalyseYield Analysis(&DataTree, true, KpiSettings::Get().GetString("PeakingBackgroundsFile"));
  Analysis.CalculateDoubleTagYields(AnalysisBinMigration.GetBinMigrationMatrix(), KpiSettings::Get().GetString("DataResultsFile"));
  Analysis.SaveDalitzDistributions(KpiSettings::Get().GetString("DalitzPlotFilename"));
  std::string FlavourTagCorrectionsFile = KpiSettings::Get().GetString("FlavourTagCorrectionsFile");
  if(FlavourTagCorrectionsFile == "None") {
    FlavourTagCorrectionsFile = "";
  }
  AnalyseFinalYields AnalysisFinalYields(KpiSettings::Get().GetString("TruthResultsFile"), KpiSettings::Get().GetString("SignalMCResultsFile"), KpiSettings::Get().GetString("DataResultsFile"), FlavourTagCorrectionsFile);
  AnalysisFinalYields.CalculateFinalYields(KpiSettings::Get().GetString("FinalYieldsFile"));
  return 0;
}
