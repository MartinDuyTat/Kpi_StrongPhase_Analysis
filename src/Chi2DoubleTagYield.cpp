// Martin Duy Tat 13th May 2021

#include<string>
#include<vector>
#include<iostream>
#include<stdexcept>
#include<iostream>
#include<numeric>
#include<utility>
#include"Chi2DoubleTagYield.h"
#include"DoubleTagMeasurement.h"
#include"TMath.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"
#include"TCanvas.h"
#include"TGraph.h"

Chi2DoubleTagYield::Chi2DoubleTagYield(bool FixNormalization, const std::string &ErrorCategory): m_FixNormalization(FixNormalization), m_ErrorCategory(ErrorCategory) {
}

void Chi2DoubleTagYield::AddMeasurement(int NBins, const std::string &K0Mode, const std::string &DataSetsToFit, const std::string &cisiHParameterFilename, const std::string &KiHParameterFilename, const std::string &DTYieldFilename) {
  if(DataSetsToFit.find(K0Mode) != std::string::npos) {
    m_Measurements.push_back(DoubleTagMeasurement(NBins, K0Mode, cisiHParameterFilename, KiHParameterFilename, DTYieldFilename));
  }
  m_cisiCovariance.AddDataset(NBins, K0Mode, cisiHParameterFilename);
}

double Chi2DoubleTagYield::operator()(const double *params) {
  double rDcosDelta = *(params + 0);
  double rDsinDelta = *(params + 1);
  std::vector<double> Normalizations(m_Measurements.size());
  for(unsigned int i = 0; i < Normalizations.size(); i++) {
    Normalizations[i] = *(params + i + 2);
  }
  double Chi2 = 0.0;
  for(unsigned int i = 0; i < m_Measurements.size(); i++) {
    Chi2 += m_Measurements[i].GetChi2(Normalizations[i], rDcosDelta, rDsinDelta, m_ErrorCategory, m_VetoBins);
  }
  return Chi2;
}

ROOT::Minuit2::Minuit2Minimizer* Chi2DoubleTagYield::RunMinimization(double &rDcosDelta, double &rDsinDelta, double &rDcosDeltaError, double &rDsinDeltaError, double &Chi2, double &Correlation) const {
  ROOT::Minuit2::Minuit2Minimizer *Minimizer = new ROOT::Minuit2::Minuit2Minimizer;
  ROOT::Math::Functor fcn(*this, 2 + m_Measurements.size());
  Minimizer->SetFunction(fcn);
  Minimizer->SetVariable(0, "rDcosDelta", 0.0, 1.0);
  Minimizer->SetVariable(1, "rDsinDelta", 0.0, 1.0);
  for(unsigned int i = 0; i < m_Measurements.size(); i++) {
    Minimizer->SetVariable(i + 2, std::string("Normalization") + std::to_string(i), 0.9, 1.1);
  }
  if(m_FixNormalization) {
    for(unsigned int i = 0; i < m_Measurements.size(); i++) {
      Minimizer->SetVariableValue(i + 2, 1.0);
      Minimizer->FixVariable(i + 2);
    }
  }
  Minimizer->Minimize();
  const double *Results = Minimizer->X();
  const double *Errors = Minimizer->Errors();
  Correlation = Minimizer->Correlation(0, 1);
  rDcosDelta = Results[0];
  rDsinDelta = Results[1];
  rDcosDeltaError = Errors[0];
  rDsinDeltaError = Errors[1];
  Chi2 = Minimizer->MinValue()/static_cast<double>(GetDegreesOfFreedom());
  return Minimizer;
}

void Chi2DoubleTagYield::MinimizeChi2(const std::string &PlotContourFilename) {
  ROOT::Minuit2::Minuit2Minimizer *Minimizer = RunMinimization(m_FittedrDcosDelta, m_FittedrDsinDelta, m_ErrorrDcosDelta, m_ErrorrDsinDelta, m_Chi2, m_Correlation);
  if(PlotContourFilename != "None") {
    DrawContours(Minimizer, PlotContourFilename);
  }
  delete Minimizer;
}

void Chi2DoubleTagYield::DrawContours(ROOT::Minuit2::Minuit2Minimizer *Minimizer, const std::string &Filename) const {
  unsigned int Npoints = 100;
  std::vector<double> x(Npoints + 1), y(Npoints + 1);
  TCanvas c("c", "c", 1200, 900);
  Minimizer->SetErrorDef(9.0);
  Minimizer->Contour(0, 1, Npoints, x.data(), y.data());
  x[Npoints] = x[0];
  y[Npoints] = y[0];
  TGraph gr1(Npoints + 1, x.data(), y.data());
  gr1.Draw("AC");
  gr1.SetTitle("r_{D}^{K#pi}cos#delta_{D}^{K#pi} vs r_{D}^{K#pi}sin#delta_{D}^{K#pi};r_{D}^{K#pi}cos#delta_{D}^{K#pi};r_{D}^{K#pi}sin#delta_{D}^{K#pi}");
  Minimizer->SetErrorDef(4.0);
  Minimizer->Contour(0, 1, Npoints, x.data(), y.data());
  x[Npoints] = x[0];
  y[Npoints] = y[0];
  TGraph gr2(Npoints + 1, x.data(), y.data());
  gr2.Draw("C");
  Minimizer->SetErrorDef(1.0);
  Minimizer->Contour(0, 1, Npoints, x.data(), y.data());
  x[Npoints] = x[0];
  y[Npoints] = y[0];
  TGraph gr3(Npoints + 1, x.data(), y.data());
  gr3.Draw("C");
  c.SaveAs(Filename.c_str());
}

void Chi2DoubleTagYield::RunSystematics(const std::string &Systematics, double &rDcosDelta_Bias, double &rDsinDelta_Bias, double &rDcosDelta_Syst, double &rDsinDelta_Syst, double &Correlation) {
  int NFits = 100000;
  if(Systematics != "Ki" && Systematics != "cisi") {
    throw std::invalid_argument("Unknown systematics");
  }
  if(Systematics != "Ki") {
    m_cisiCovariance.PrepareCholesky();
  }
  std::cout << "Running many fits for " << Systematics << " systematics\n";
  std::vector<double> rDcosDeltaFit(NFits), rDsinDeltaFit(NFits);
  for(int i = 0; i < NFits; i++) {
    if(i%(NFits/10) == 0) {
      std::cout << "Fit number " << i << "\n";
    }
    if(Systematics != "Ki") {
      m_cisiCovariance.Smear();
    }
    for(auto& Measurement : m_Measurements) {
      if(Systematics == "Ki") {
	Measurement.SmearKi();
      } else {
	std::string Mode = Measurement.Mode();	
	Measurement.Smearcisi(m_cisiCovariance.GetSmearing(Mode));
      }
    }
    double dummy;
    RunMinimization(rDcosDeltaFit[i], rDsinDeltaFit[i], dummy, dummy, dummy, dummy);
  }
  for(auto& Measurement : m_Measurements) {
    Measurement.RemoveSmearing();
  }
  rDcosDelta_Bias = TMath::Mean(rDcosDeltaFit.begin(), rDcosDeltaFit.end());
  rDsinDelta_Bias = TMath::Mean(rDsinDeltaFit.begin(), rDsinDeltaFit.end());
  rDcosDelta_Syst = TMath::RMS(rDcosDeltaFit.begin(), rDcosDeltaFit.end());
  rDsinDelta_Syst = TMath::RMS(rDsinDeltaFit.begin(), rDsinDeltaFit.end());
  std::transform(rDcosDeltaFit.begin(), rDcosDeltaFit.end(), rDcosDeltaFit.begin(), [&rDcosDelta_Bias] (double a) { return a - rDcosDelta_Bias; });
  std::transform(rDsinDeltaFit.begin(), rDsinDeltaFit.end(), rDsinDeltaFit.begin(), [&rDsinDelta_Bias] (double a) { return a - rDsinDelta_Bias; });
  Correlation = std::inner_product(rDcosDeltaFit.begin(), rDcosDeltaFit.end(), rDsinDeltaFit.begin(), 0.0)/((NFits - 1)*rDcosDelta_Syst*rDsinDelta_Syst);
}

double Chi2DoubleTagYield::GetFittedrDcosDelta() const {
  return m_FittedrDcosDelta;
}

double Chi2DoubleTagYield::GetErrorrDcosDelta() const {
  return m_ErrorrDcosDelta;
}

double Chi2DoubleTagYield::GetFittedrDsinDelta() const {
  return m_FittedrDsinDelta;
}

double Chi2DoubleTagYield::GetErrorrDsinDelta() const {
  return m_ErrorrDsinDelta;
}

double Chi2DoubleTagYield::GetChi2PerDegreesOfFreedom() const {
  return m_Chi2;
}

int Chi2DoubleTagYield::GetDegreesOfFreedom() const {
  int DOF = 0;
  for(auto iter = m_Measurements.begin(); iter != m_Measurements.end(); iter++) {
    DOF += 2*iter->GetNBins();
  }
  if(!m_FixNormalization) {
    DOF -= m_Measurements.size();
  }
  return DOF - 2;
}

double Chi2DoubleTagYield::GetCorrelation() const {
  return m_Correlation;
}

void Chi2DoubleTagYield::SetVetoBins(const std::vector<std::pair<std::string, int>> &VetoBins) {
  m_VetoBins = VetoBins;
}
