// Martin Duy Tat 19th May 2021

#include<stdexcept>
#include<utility>
#include<string>
#include"BinningScheme.h"
#include"TFile.h"
#include"TStyle.h"

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

int BinningScheme::GetMappedBinNumber(double M2Plus, double M2Minus, int KCharge) const {
  int CP = M2Plus > M2Minus ? +1 : -1;
  // First get initial bin number
  int x = m_BinningScheme->GetXaxis()->FindBin(M2Plus);
  int y = m_BinningScheme->GetYaxis()->FindBin(M2Minus);
  int i = 1;
  // Loop in circles around this point until reaching the Dalitz boundary
  while(true) {
    if(i > 1000) {
      throw std::runtime_error("Dalitz point too far outside phase space");
    }
    for(int j = 0; j <= i; j++) {
      // All possible combinations of displacements at the same distance
      std::vector<std::pair<int, int>> BinList{{i, j}, {i, -j}, {-i, j}, {-i, -j}, {j, i}, {j, -i}, {-j, i}, {-j, -i}};
      for(auto iter = BinList.begin(); iter != BinList.end(); iter++) {
        int NewBin = static_cast<int>(m_BinningScheme->GetBinContent(x + iter->first, y + iter->second));
	// Once we reach the Dalitz boundary the bin number is non-zero
        if(NewBin != 0) {
	  NewBin *= KCharge*CP;
	  return NewBin;
        }
      }
    }
    i++;
  }
}

int BinningScheme::GetNumberBins() const {
  return m_NumberBins;
}

void BinningScheme::Draw(const std::string &DrawOptions, const std::string &PlotTitle) const {
  int palette[3] = {kWhite, kOrange - 2, kOrange - 3};
  gStyle->SetPalette(3, palette);
  m_BinningScheme->SetTitle(PlotTitle.c_str());
  m_BinningScheme->Draw(DrawOptions.c_str());
}
