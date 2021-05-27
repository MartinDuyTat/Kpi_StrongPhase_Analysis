// Martin Duy Tat 19th May 2021
/**
 * AnalyseBackgrounds is an application for determining the binned yields of various backgrounds after running through TopoAna
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
#include"AnalysePeakingBackgrounds.h"

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 1 input argument\n";
    return 0;
  }
  KpiSettings::Get().Initialize(std::string(argv[1]));
  TreeWrapper TopoAnaTree(KpiSettings::Get().GetString("TopoAnaFiles"), KpiSettings::Get().GetString("TreeName"), KpiSettings::Get().GetString("DataCutFile"), "TopoAna");
  AnalysePeakingBackgrounds AnalysisBackgrounds(&TopoAnaTree, KpiSettings::Get().GetString("BackgroundTopologiesFile"));
  AnalysisBackgrounds.CalculatePeakingBackgrounds(KpiSettings::Get().GetString("PeakingBackgroundsFile"), KpiSettings::Get().GetDouble("LuminosityScale"));
  return 0;
}
