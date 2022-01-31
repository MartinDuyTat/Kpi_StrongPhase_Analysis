// Martin Duy Tat 21st July 2021

#include<string>
#include"PresentYields.h"
#include"BinVector.h"
#include"Bes3plotstyle.h"
#include"TStyle.h"
#include"TCanvas.h"
#include"TH1D.h"
#include"TLegend.h"
#include"TLine.h"
#include"TLatex.h"
#include"TPad.h"

PresentYields::PresentYields(const std::string &K0pipiKpiYieldFile, const std::string &K0KKKpiYieldFile, const std::string &K0pipicisiHadronicFile, const std::string &K0pipiKiHadronicFile, const std::string &K0KKcisiHadronicFile, const std::string &K0KKKiHadronicFile, const std::string &K0Mode, double rDcosDelta, double rDsinDelta, bool DrawK0KK, bool PlotPulls): m_K0pipiBins(8), m_K0KKBins(2), m_K0pipiDoubleTagYield(K0pipiKpiYieldFile, m_K0pipiBins), m_K0KKDoubleTagYield(K0KKKpiYieldFile, m_K0KKBins), m_K0pipiHadronicParameters(K0pipicisiHadronicFile, K0pipiKiHadronicFile, m_K0pipiBins, K0Mode), m_K0KKHadronicParameters(K0KKcisiHadronicFile, K0KKKiHadronicFile, m_K0KKBins, K0Mode), m_K0Mode(K0Mode), m_rDcosDelta(rDcosDelta), m_rDsinDelta(rDsinDelta), m_PlotPulls(PlotPulls) {
  if(!DrawK0KK) {
    m_K0KKBins = 0;
  }
}

void PresentYields::PlotYieldPresentation(const std::string &Filename) const {
  //gStyle->SetErrorX(0);
  //gStyle->SetOptStat(0);
  SetStyle();
  SetPrelimStyle();
  TCanvas c("c", "c", 1500, 1000);
  TPad Pad1("p1", "p1", 0.0, m_PlotPulls ? 0.2 : 0.0, 1.0, 1.0), Pad2("p2", "p2", 0.0, 0.0, 1.0, m_PlotPulls ? 0.2 : 0.0);
  //Pad1.SetBottomMargin(0.1);
  //Pad2.SetTopMargin(0);
  //Pad2.SetBottomMargin(0.1);
  Pad1.Draw();
  Pad2.Draw();
  TH1D h1("h1", "Fit", 2*m_K0pipiBins + 2*m_K0KKBins, 0, 2*m_K0pipiBins + 2*m_K0KKBins);
  TH1D h2("h2", "Uncorrelated prediction", 2*m_K0pipiBins + 2*m_K0KKBins, 0, 2*m_K0pipiBins + 2*m_K0KKBins);
  TH1D h3("h3", "Measurement", 2*m_K0pipiBins + 2*m_K0KKBins, 0, 2*m_K0pipiBins + 2*m_K0KKBins);
  TH1D h4("h4", "", 2*m_K0pipiBins + 2*m_K0KKBins, 0, 2*m_K0pipiBins + 2*m_K0KKBins);
  h1.SetLineWidth(3);
  h2.SetLineWidth(3);
  //h2.SetLineStyle(kDashed);
  h3.SetLineWidth(3);
  h3.SetLineColor(kBlack);
  h4.SetLineColor(kBlack);
  h4.SetLineWidth(3);
  BinVector<double> PredictedK0pipiYield, PredictedK0KKYield, dummy;
  m_K0pipiHadronicParameters.CalculateNormalizedYields(1, m_rDcosDelta, m_rDsinDelta, PredictedK0pipiYield, dummy, dummy);
  m_K0KKHadronicParameters.CalculateNormalizedYields(1, m_rDcosDelta, m_rDsinDelta, PredictedK0KKYield, dummy, dummy);
  for(int i = 1; i <= m_K0pipiBins; i++) {
    h1.SetBinContent(i, PredictedK0pipiYield[i]);
    h2.SetBinContent(i, m_K0pipiHadronicParameters.GetKi(i));
    h3.SetBinContent(i, m_K0pipiDoubleTagYield.GetYield(i));
    h3.SetBinError(i, m_K0pipiDoubleTagYield.GetYieldError(i));
    h4.SetBinContent(i, (m_K0pipiDoubleTagYield.GetYield(i) - PredictedK0pipiYield[i])/m_K0pipiDoubleTagYield.GetYieldError(i));
    h4.SetBinError(i, 1.0);
    h1.GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
    h4.GetXaxis()->SetBinLabel(i, std::to_string(i).c_str());
  }
  for(int i = 1; i <= m_K0pipiBins; i++) {
    h1.SetBinContent(i + m_K0pipiBins, PredictedK0pipiYield[-i]);
    h2.SetBinContent(i + m_K0pipiBins, m_K0pipiHadronicParameters.GetKi(-i));
    h3.SetBinContent(i + m_K0pipiBins, m_K0pipiDoubleTagYield.GetYield(-i));
    h3.SetBinError(i + m_K0pipiBins, m_K0pipiDoubleTagYield.GetYieldError(-i));
    h4.SetBinContent(i + m_K0pipiBins, (m_K0pipiDoubleTagYield.GetYield(-i) - PredictedK0pipiYield[-i])/m_K0pipiDoubleTagYield.GetYieldError(-i));
    h4.SetBinError(i + m_K0pipiBins, 1.0);
    h1.GetXaxis()->SetBinLabel(i + m_K0pipiBins, std::to_string(-i).c_str());
    h4.GetXaxis()->SetBinLabel(i + m_K0pipiBins, std::to_string(-i).c_str());
  }
  for(int i = 1; i <= m_K0KKBins; i++) {
    h1.SetBinContent(i + 2*m_K0pipiBins, PredictedK0KKYield[i]/2);
    h2.SetBinContent(i + 2*m_K0pipiBins, m_K0KKHadronicParameters.GetKi(i)/2);
    h3.SetBinContent(i + 2*m_K0pipiBins, m_K0KKDoubleTagYield.GetYield(i)/2);
    h3.SetBinError(i + 2*m_K0pipiBins, m_K0KKDoubleTagYield.GetYieldError(i)/2);
    h4.SetBinContent(i + 2*m_K0pipiBins, (m_K0KKDoubleTagYield.GetYield(i) - PredictedK0KKYield[i])/m_K0KKDoubleTagYield.GetYieldError(i));
    h4.SetBinError(i + 2*m_K0pipiBins, 1.0);
    h1.GetXaxis()->SetBinLabel(i + 2*m_K0pipiBins, std::to_string(i).c_str());
    h4.GetXaxis()->SetBinLabel(i + 2*m_K0pipiBins, std::to_string(i).c_str());
  }
  for(int i = 1; i <= m_K0KKBins; i++) {
    h1.SetBinContent(i + 2*m_K0pipiBins + m_K0KKBins, PredictedK0KKYield[-i]/2);
    h2.SetBinContent(i + 2*m_K0pipiBins + m_K0KKBins, m_K0KKHadronicParameters.GetKi(-i)/2);
    h3.SetBinContent(i + 2*m_K0pipiBins + m_K0KKBins, m_K0KKDoubleTagYield.GetYield(-i)/2);
    h3.SetBinError(i + 2*m_K0pipiBins + m_K0KKBins, m_K0KKDoubleTagYield.GetYieldError(-i)/2);
    h4.SetBinContent(i + 2*m_K0pipiBins + m_K0KKBins, (m_K0KKDoubleTagYield.GetYield(-i) - PredictedK0KKYield[-i])/m_K0KKDoubleTagYield.GetYieldError(-i));
    h4.SetBinError(i + 2*m_K0pipiBins + m_K0KKBins, 1.0);
    h1.GetXaxis()->SetBinLabel(i + 2*m_K0pipiBins + m_K0KKBins, std::to_string(-i).c_str());
    h4.GetXaxis()->SetBinLabel(i + 2*m_K0pipiBins + m_K0KKBins, std::to_string(-i).c_str());
  }
  Pad1.cd();
  TLegend legend(0.57,0.7,0.95,0.90);
  legend.AddEntry(&h3, "", "LEP");
  legend.AddEntry(&h1, "", "L");
  legend.AddEntry(&h2, "", "L");
  legend.SetFillColor(0);
  legend.SetTextFont(42);
  std::string K0ModeText = m_K0Mode == "KS" ? std::string("K_{S}") : std::string("K_{L}");
  //std::string Title = K0ModeText + std::string("hh vs K#pi double tag yield measurement and prediction;Bin;Normalised yield");
  //h1.SetTitle(Title.c_str());
  NameAxes(&h1, const_cast<char*>(std::string("Bin number").c_str()), const_cast<char*>(std::string("Normalised yield").c_str()));
  FormatData(&h3);
  if(m_K0Mode == "KS" && !m_PlotPulls) {
    h1.GetXaxis()->SetLabelOffset(999);
    h1.GetXaxis()->SetLabelSize(0);
    h1.SetXTitle("");
  } else if(m_K0Mode == "KL" && !m_PlotPulls) {
    h1.SetYTitle("");
  }
  h1.SetLineColor(2);
  h1.SetLineStyle(1);
  h2.SetLineColor(4);
  h2.SetLineStyle(2);
  h1.Draw("HIST");
  h2.Draw("HIST SAME");
  h3.Draw("P E1 SAME");
  TLine Line(2*m_K0pipiBins, 0.0, 2*m_K0pipiBins, h1.GetBinContent(h1.GetMaximumBin()));
  Line.SetLineStyle(kDashed);
  Line.SetLineColor(kBlack);
  if(m_K0KKBins != 0) {
    Line.Draw("SAME");
  }
  if(m_K0Mode == "KS") {
    legend.Draw();
  }
  TLatex K0pipiText;
  K0pipiText.SetTextSize(0.069);
  K0pipiText.SetNDC();
  K0pipiText.DrawLatex(0.67, 0.6, ("#it{#bf{D #rightarrow " + K0ModeText + std::string("^{0}#pi^{+}#pi^{#minus}}}")).c_str());
  TLatex K0KKText;
  K0KKText.SetTextSize(0.069);
  K0KKText.SetNDC();
  if(m_K0KKBins != 0) {
    K0KKText.DrawLatex(0.76, 0.6, ("#it{#bf{D #rightarrow " + K0ModeText + std::string("^{0}K^{+}K^{#minus}}}")).c_str());
  }
  if(m_PlotPulls) {
    Pad2.cd();
    h4.GetXaxis()->SetLabelFont(0);
    h4.GetXaxis()->SetLabelSize(0);
    h4.GetYaxis()->SetLabelFont(62);
    h4.GetYaxis()->SetLabelSize(0.1);
    h4.SetMaximum(4.0);
    h4.SetMinimum(-4.0);
    h4.Draw("P E1");
    TLine Line2(0.0, 0.0, 2*m_K0pipiBins + 2*m_K0KKBins, 0.0);
    Line2.SetLineStyle(kDashed);
    Line2.SetLineColor(kBlack);
    Line2.Draw("SAME");
  }
  c.SaveAs(Filename.c_str());
}
