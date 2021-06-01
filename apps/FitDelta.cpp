// Martin Duy Tat 18th May 2021
/**
 * This application fits the measurements of double tag yields to obtain \f$r_D^{K\pi}\cos(\delta_D^{K\pi})\f$ and \f$r_D^{K\pi}\sin(\delta_D^{K\pi})\f$
 * @param 1 Filename of settings file
 */

#include<iostream>
#include<string>
#include"Settings.h"
#include"DeltaKpiFitter.h"

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 1 input argument\n";
    return 0;
  }
  KpiSettings::Get().Initialize(std::string(argv[1]));
  DeltaKpiFitter Fitter(KpiSettings::Get().GetString("Datasets"));
  Fitter.RunFit(KpiSettings::Get().GetString("ResultsFile"));
  return 0;
}