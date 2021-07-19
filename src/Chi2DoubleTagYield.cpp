// Martin Duy Tat 13th May 2021

#include"Chi2DoubleTagYield.h"
#include"DoubleTagMeasurement.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"

Chi2DoubleTagYield::Chi2DoubleTagYield(bool FixNormalization): m_FixNormalization(FixNormalization) {
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
    Chi2 += m_Measurements[i].GetChi2(Normalizations[i], rDcosDelta, rDsinDelta);
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
