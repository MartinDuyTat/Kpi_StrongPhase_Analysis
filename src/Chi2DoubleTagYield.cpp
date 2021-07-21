// Martin Duy Tat 13th May 2021

#include<string>
#include<vector>
#include<iostream>
#include"Chi2DoubleTagYield.h"
#include"DoubleTagMeasurement.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"
#include"TCanvas.h"
#include"TGraph.h"

Chi2DoubleTagYield::Chi2DoubleTagYield(bool FixNormalization, const std::string &ErrorCategory): m_FixNormalization(FixNormalization), m_ErrorCategory(ErrorCategory) {
}

void Chi2DoubleTagYield::AddMeasurement(const DoubleTagMeasurement &Measurement) {
  m_Measurements.push_back(Measurement);
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
    Chi2 += m_Measurements[i].GetChi2(Normalizations[i], rDcosDelta, rDsinDelta, m_ErrorCategory);
  }
  return Chi2;
}

void Chi2DoubleTagYield::MinimizeChi2() {
  ROOT::Minuit2::Minuit2Minimizer Minimizer;
  ROOT::Math::Functor fcn(*this, 2 + m_Measurements.size());
  Minimizer.SetFunction(fcn);
  Minimizer.SetVariable(0, "rDcosDelta", 0.0, 1.0);
  Minimizer.SetVariable(1, "rDsinDelta", 0.0, 1.0);
  for(unsigned int i = 0; i < m_Measurements.size(); i++) {
    Minimizer.SetVariable(i + 2, std::string("Normalization") + std::to_string(i), 0.9, 1.1);
  }
  if(m_FixNormalization) {
    for(unsigned int i = 0; i < m_Measurements.size(); i++) {
      Minimizer.SetVariableValue(i + 2, 1.0);
      Minimizer.FixVariable(i + 2);
    }
  }
  Minimizer.Minimize();
  const double *Results = Minimizer.X();
  const double *Errors = Minimizer.Errors();
  m_FittedrDcosDelta = Results[0];
  m_FittedrDsinDelta = Results[1];
  m_ErrorrDcosDelta = Errors[0];
  m_ErrorrDsinDelta = Errors[1];
  m_Chi2 = Minimizer.MinValue()/static_cast<double>(GetDegreesOfFreedom());
  std::cout << "Plot contours? ";
  std::string Answer;
  std::cin >> Answer;
  if(Answer == "yes") {
    std::cout << "Filename: ";
    std::string Filename;
    std::cin >> Filename;
    DrawContours(&Minimizer, Filename);
  }
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
