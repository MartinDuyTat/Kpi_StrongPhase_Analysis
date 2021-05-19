// Martin Duy Tat 19th May 2021

#include"BinningScheme.h"
#include"TFile.h"

BinningScheme::BinningScheme(): m_NumberBins(2) {
  TFile BinningFile((std::string(BINNING_SCHEME_DIR) + std::string("KsKK_2bins.root")).c_str(), "READ");
  BinningFile.GetObject("dkpp_bin_h", m_BinningScheme);
  m_BinningScheme->SetDirectory(0);
  BinningFile.Close();
}

int BinningScheme::GetBinNumber(double M2Plus, double M2Minus, int KCharge) const {
  Float_t BinNumberFloat = m_BinningScheme->GetBinContent(m_BinningScheme->GetXaxis()->FindBin(M2Plus), m_BinningScheme->GetYaxis()->FindBin(M2Minus));
  int CP = M2Plus > M2Minus ? +1 : -1;
  int BinNumber = static_cast<int>(BinNumberFloat);
  BinNumber *= KCharge*CP;
  return BinNumber;
}

int BinningScheme::GetNumberBins() const {
  return m_NumberBins;
}
