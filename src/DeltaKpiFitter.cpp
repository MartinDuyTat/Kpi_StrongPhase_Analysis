// Martin Duy Tat 18th May 2021

#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include"DeltaKpiFitter.h"

DeltaKpiFitter::DeltaKpiFitter(const std::string &Filename, const std::string &DataSetsToFit, bool FixNormalization, const std::string &ErrorCategory): m_Measurements(Chi2DoubleTagYield(FixNormalization, ErrorCategory)) {
  std::ifstream Infile(Filename);
  std::string line;
  int N = 0;
  while(std::getline(Infile, line)) {
    if(line[0] == '#') {
      continue;
    }
    std::stringstream ss(line);
    int NBins;
    std::string K0Mode, HParameterFilename, DTYieldFilename;
    ss >> NBins >> K0Mode >> HParameterFilename >> DTYieldFilename;
    m_Measurements.AddMeasurement(NBins, K0Mode, DataSetsToFit, HParameterFilename, DTYieldFilename);
    N++;
  }
  Infile.close();
  if(N == 0) {
    std::cout << "No datasets inputted\n";
  } else if(N == 1) {
    std::cout << "1 dataset added\n";
  } else {
    std::cout << N << " datasets added\n";
  }
}

void DeltaKpiFitter::RunFit(const std::string &Filename, const std::string &PlotContourFilename, bool RunKiSystematics, bool RuncisiSystematics) {
  m_Measurements.MinimizeChi2(PlotContourFilename);
  std::ofstream Outfile(Filename);
  Outfile << "Statistical fit\n";
  Outfile << "Chi2/DOF: " << m_Measurements.GetChi2PerDegreesOfFreedom() << "\n";
  Outfile << "r_D*cos(delta_D): " << m_Measurements.GetFittedrDcosDelta() << " \u00B1 " << m_Measurements.GetErrorrDcosDelta() << "\n";
  Outfile << "r_D*sin(delta_D): " << m_Measurements.GetFittedrDsinDelta() << " \u00B1 " << m_Measurements.GetErrorrDsinDelta() << "\n";
  Outfile << "Correlation: " << m_Measurements.GetCorrelation() << "\n";
  if(RunKiSystematics) {
    double Ki_rDcosDelta_Bias, Ki_rDsinDelta_Bias, Ki_rDcosDelta_Syst, Ki_rDsinDelta_Syst, Ki_Correlation;
    m_Measurements.RunSystematics("Ki", Ki_rDcosDelta_Bias, Ki_rDsinDelta_Bias, Ki_rDcosDelta_Syst, Ki_rDsinDelta_Syst, Ki_Correlation);
    Outfile << "Ki systematics\n";
    Outfile << "r_D*cos(delta_D): " << Ki_rDcosDelta_Bias << " \u00B1 " << Ki_rDcosDelta_Syst << "\n";
    Outfile << "r_D*sin(delta_D): " << Ki_rDsinDelta_Bias << " \u00B1 " << Ki_rDsinDelta_Syst << "\n";
    Outfile << "Correlation: " << Ki_Correlation << "\n";
  }
  if(RuncisiSystematics) {
    double cisi_rDcosDelta_Bias, cisi_rDsinDelta_Bias, cisi_rDcosDelta_Syst, cisi_rDsinDelta_Syst, cisi_Correlation;
    Outfile << "cisi systematics\n";
    m_Measurements.RunSystematics("cisi", cisi_rDcosDelta_Bias, cisi_rDsinDelta_Bias, cisi_rDcosDelta_Syst, cisi_rDsinDelta_Syst, cisi_Correlation);
    Outfile << cisi_rDcosDelta_Bias << " \u00B1 " << cisi_rDcosDelta_Syst << "\n";
    Outfile << cisi_rDsinDelta_Bias << " \u00B1 " << cisi_rDsinDelta_Syst << "\n";
    Outfile << "Correlation: " << cisi_Correlation << "\n";
  }
  Outfile.close();
}
