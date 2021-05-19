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

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 1 input argument\n";
    return 0;
  }
  KpiSettings::Get().Initialize(std::string(argv[1]));
  TreeWrapper Tree(KpiSettings::Get().GetString("DataFiles"), KpiSettings::Get().GetString("TreeName"), KpiSettings::Get().GetString("CutFile"), "Data");
  AnalyseYield Analysis(&Tree);
  Analysis.CalculateDoubleTagYields(KpiSettings::Get().GetString("ResultsFile"));
  return 0;
}
