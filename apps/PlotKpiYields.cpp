// Martin Duy Tat 21st July 2021
/**
 * PlotKpiYields is an application that plots the Kpi double tag yields and the predicted for each bin
 * @param 1 Filename of settings file
 */

#include<iostream>
#include"PresentYields.h"
#include"Settings.h"
#include"TMath.h"

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 1 input argument\n";
    return 0;
  }
  std::cout << "Plotting Kpi bin yields and prediction\n";
  KpiSettings::Get().Initialize(std::string(argv[1]));
  std::cout << "Loading data...\n";
  PresentYields PresentYield(KpiSettings::Get().GetString("K0pipiKpiYieldFile"), KpiSettings::Get().GetString("K0KKKpiYieldFile"), KpiSettings::Get().GetString("K0pipicisiHadronicFile"), KpiSettings::Get().GetString("K0pipiKiHadronicFile"), KpiSettings::Get().GetString("K0KKcisiHadronicFile"), KpiSettings::Get().GetString("K0KKKiHadronicFile"), KpiSettings::Get().GetString("K0Mode"), -0.0562213, -0.011149, KpiSettings::Get().GetBool("DrawK0KK"));
  std::cout << "Data ready\n";
  std::cout << "Plotting...\n";
  PresentYield.PlotYieldPresentation(KpiSettings::Get().GetString("PlotName"));
  std::cout << "Plot saved\n";
  return 0;
}
