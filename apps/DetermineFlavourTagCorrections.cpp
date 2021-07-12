// Martin Duy Tat 9th July 2021
/**
 * DetermineFlavourTagCorrections is an application that calculates the flavour tag corrections due to DCS decays
 * @param 1 Filename of settings file
 */

#include<iostream>
#include"FlavourCorrection.h"
#include"Settings.h"

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 1 input argument\n";
    return 0;
  }
  KpiSettings::Get().Initialize(std::string(argv[1]));
  std::cout << "Calculating flavour tag corrections...\n";
  FlavourCorrection FlavourCorrections(KpiSettings::Get().GetInt("Bins"), KpiSettings::Get().GetString("Mode"), KpiSettings::Get().GetString("StrongPhasesFromModelFile"), KpiSettings::Get().GetString("HadronicParametersFile"));
  FlavourCorrections.SaveCorrections(KpiSettings::Get().GetString("FlavourTagCorrectionsFile"));
  std::cout << "Flavour tag corrections ready!\n";
  return 0;
}
