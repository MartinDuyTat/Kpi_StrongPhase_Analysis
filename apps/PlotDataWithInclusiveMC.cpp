// Martin Duy Tat 29th May 2021
/**
 * PlotDataWithInclusiveMC is an application that plots the inclusive MC samples and compares it with the data
 * There is also an option to remove the datapoints, or to remove certain peaking backgrounds
 * It also counts the number of events in each region, which can be used to calculate the flat backgrounds in each region after all peaking backgrounds are removed
 * @param 1 Filename of settings file
 */

#include<iostream>
#include<fstream>
#include<map>
#include"Settings.h"
#include"TFile.h"
#include"TTree.h"
#include"TChain.h"
#include"TH1D.h"
#include"TCut.h"
#include"THStack.h"
#include"TCanvas.h"
#include"TPad.h"
#include"TStyle.h"

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 1 input argument\n";
    return 0;
  }
  KpiSettings::Get().Initialize(std::string(argv[1]));
  std::ifstream CutFile(KpiSettings::Get().GetString("DataCutFile"));
  TCut Cuts;
  std::string line;
  while(std::getline(CutFile, line)) {
    Cuts = Cuts && TCut(line.c_str());
  }
  CutFile.close();
  std::string PlotVariable = KpiSettings::Get().GetString("PlotVariable");
  std::string LowSidebandCut = KpiSettings::Get().GetString("LowSidebandCut");
  std::string SignalRegionCut = KpiSettings::Get().GetString("SignalRegionCut");
  std::string UpperSidebandCut = KpiSettings::Get().GetString("UpperSidebandCut");
  std::map<char, double> Yields{{'L', 0.0}, {'S', 0.0}, {'H', 0.0}};
  // qq2010
  TFile f1(KpiSettings::Get().GetString("qq2010File").c_str(), "READ");
  TTree *t1 = nullptr;
  f1.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t1);
  std::string Cut1 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("qq2010LuminosityScale");
  TH1D h1("h1", "q#bar{q}", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t1->Draw((PlotVariable + std::string(" >> h1")).c_str(), Cut1.c_str(), "goff");
  Yields['L'] += t1->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("qq2010LuminosityScale");
  Yields['S'] += t1->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("qq2010LuminosityScale");
  Yields['H'] += t1->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("qq2010LuminosityScale");
  // qq2011
  TFile f2(KpiSettings::Get().GetString("qq2011File").c_str(), "READ");
  TTree *t2 = nullptr;
  f2.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t2);
  std::string Cut2 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("qq2011LuminosityScale");
  TH1D h2("h2", "q#bar{q}", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t2->Draw((PlotVariable + std::string(" >> h2")).c_str(), Cut2.c_str(), "goff");
  Yields['L'] += t2->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("qq2011LuminosityScale");
  Yields['S'] += t2->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("qq2011LuminosityScale");
  Yields['H'] += t2->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("qq2011LuminosityScale");
  h1.Add(&h2);
  // DpDm2010
  TFile f3(KpiSettings::Get().GetString("DpDm2010File").c_str(), "READ");
  TTree *t3 = nullptr;
  f3.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t3);
  std::string Cut3 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("DpDm2010LuminosityScale");
  TH1D h3("h3", "D^{+}D^{#minus}", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t3->Draw((PlotVariable + std::string(" >> h3")).c_str(), Cut3.c_str(), "goff");
  Yields['L'] += t3->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("DpDm2010LuminosityScale");
  Yields['S'] += t3->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("DpDm2010LuminosityScale");
  Yields['H'] += t3->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("DpDm2010LuminosityScale");
  // DpDm2011
  TFile f4(KpiSettings::Get().GetString("DpDm2011File").c_str(), "READ");
  TTree *t4 = nullptr;
  f4.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t4);
  std::string Cut4 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("DpDm2011LuminosityScale");
  TH1D h4("h4", "D^{+}D^{#minus}", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t4->Draw((PlotVariable + std::string(" >> h4")).c_str(), Cut4.c_str(), "goff");
  Yields['L'] += t4->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("DpDm2011LuminosityScale");
  Yields['S'] += t4->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("DpDm2011LuminosityScale");
  Yields['H'] += t4->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("DpDm2011LuminosityScale");
  h3.Add(&h4);
  // tautau2010
  TFile f5(KpiSettings::Get().GetString("tautau2010File").c_str(), "READ");
  TTree *t5 = nullptr;
  f5.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t5);
  std::string Cut5 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("tautau2010LuminosityScale");
  TH1D h5("h5", "#tau#tau", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t5->Draw((PlotVariable + std::string(" >> h5")).c_str(), Cut5.c_str(), "goff");
  Yields['L'] += t5->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("tautau2010LuminosityScale");
  Yields['S'] += t5->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("tautau2010LuminosityScale");
  Yields['H'] += t5->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("tautau2010LuminosityScale");
  // tautau2011
  TFile f6(KpiSettings::Get().GetString("tautau2011File").c_str(), "READ");
  TTree *t6 = nullptr;
  f6.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t6);
  std::string Cut6 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("tautau2011LuminosityScale");
  TH1D h6("h6", "#tau#tau", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t6->Draw((PlotVariable + std::string(" >> h6")).c_str(), Cut6.c_str(), "goff");
  Yields['L'] += t6->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("tautau2011LuminosityScale");
  Yields['S'] += t6->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("tautau2011LuminosityScale");
  Yields['H'] += t6->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("tautau2011LuminosityScale");
  h5.Add(&h6);
  // RR2S2010
  TFile f7(KpiSettings::Get().GetString("RR2S2010File").c_str(), "READ");
  TTree *t7 = nullptr;
  f7.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t7);
  std::string Cut7 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("RR2S2010LuminosityScale");
  TH1D h7("h7", "#psi''", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t7->Draw((PlotVariable + std::string(" >> h7")).c_str(), Cut7.c_str(), "goff");
  Yields['L'] += t7->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("RR2S2010LuminosityScale");
  Yields['S'] += t7->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("RR2S2010LuminosityScale");
  Yields['H'] += t7->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("RR2S2010LuminosityScale");
  // RR2S2011
  TFile f8(KpiSettings::Get().GetString("RR2S2011File").c_str(), "READ");
  TTree *t8 = nullptr;
  f8.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t8);
  std::string Cut8 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("RR2S2011LuminosityScale");
  TH1D h8("h8", "#psi''", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t8->Draw((PlotVariable + std::string(" >> h8")).c_str(), Cut8.c_str(), "goff");
  Yields['L'] += t8->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("RR2S2011LuminosityScale");
  Yields['S'] += t8->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("RR2S2011LuminosityScale");
  Yields['H'] += t8->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("RR2S2011LuminosityScale");
  h7.Add(&h8);
  // RR1S2010
  TFile f9(KpiSettings::Get().GetString("RR1S2010File").c_str(), "READ");
  TTree *t9 = nullptr;
  f9.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t9);
  std::string Cut9 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("RR1S2010LuminosityScale");
  TH1D h9("h9", "#psi'", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t9->Draw((PlotVariable + std::string(" >> h9")).c_str(), Cut9.c_str(), "goff");
  Yields['L'] += t9->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("RR1S2010LuminosityScale");
  Yields['S'] += t9->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("RR1S2010LuminosityScale");
  Yields['H'] += t9->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("RR1S2010LuminosityScale");
  // RR1S2011
  TFile f10(KpiSettings::Get().GetString("RR1S2011File").c_str(), "READ");
  TTree *t10 = nullptr;
  f10.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t10);
  std::string Cut10 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("RR1S2011LuminosityScale");
  TH1D h10("h10", "#psi'", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t10->Draw((PlotVariable + std::string(" >> h10")).c_str(), Cut10.c_str(), "goff");
  Yields['L'] += t10->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("RR1S2011LuminosityScale");
  Yields['S'] += t10->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("RR1S2011LuminosityScale");
  Yields['H'] += t10->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("RR1S2011LuminosityScale");
  h9.Add(&h10);
  // nonDD2010
  TFile f11(KpiSettings::Get().GetString("nonDD2010File").c_str(), "READ");
  TTree *t11 = nullptr;
  f11.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t11);
  std::string Cut11 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("nonDD2010LuminosityScale");
  TH1D h11("h11", "non D#bar{D}", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t11->Draw((PlotVariable + std::string(" >> h11")).c_str(), Cut11.c_str(), "goff");
  Yields['L'] += t11->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("nonDD2010LuminosityScale");
  Yields['S'] += t11->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("nonDD2010LuminosityScale");
  Yields['H'] += t11->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("nonDD2010LuminosityScale");
  // nonDD2011
  TFile f12(KpiSettings::Get().GetString("nonDD2011File").c_str(), "READ");
  TTree *t12 = nullptr;
  f12.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t12);
  std::string Cut12 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("nonDD2011LuminosityScale");
  TH1D h12("h12", "non D#bar{D}", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t12->Draw((PlotVariable + std::string(" >> h12")).c_str(), Cut12.c_str(), "goff");
  Yields['L'] += t12->GetEntries(Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("nonDD2011LuminosityScale");
  Yields['S'] += t12->GetEntries(Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("nonDD2011LuminosityScale");
  Yields['H'] += t12->GetEntries(Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("nonDD2011LuminosityScale");
  h11.Add(&h12);
  // D0D0
  TFile f13(KpiSettings::Get().GetString("D0D0File").c_str(), "READ");
  TTree *t13 = nullptr;
  f13.GetObject(KpiSettings::Get().GetString("TreeName").c_str(), t13);
  std::string Cut13 = std::string("(") + std::string(Cuts.GetTitle()) + std::string(")*1.0/") + KpiSettings::Get().GetString("D0D0LuminosityScale");
  TCut PeakingCut;
  if(KpiSettings::Get().GetString("RemovePeaking") == "True") {
    Cut13 += std::string("*(");
    std::ifstream iDcyTrFile(KpiSettings::Get().GetString("iDcyTrNumbersFile"));
    int iDcyTr;
    while(iDcyTrFile >> iDcyTr) {
      PeakingCut = PeakingCut && TCut((std::string("iDcyTr != ") + std::to_string(iDcyTr)).c_str());
      Cut13 += std::string("iDcyTr != ") + std::to_string(iDcyTr) + std::string(" &&");
    }
    Cut13 += std::string(" iDcyTr >= 0)");
    iDcyTrFile.close();
  }
  TH1D h13("h13", "D^{0}#bar{D^{0}}", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  t13->Draw((PlotVariable + std::string(" >> h13")).c_str(), Cut13.c_str(), "goff");
  Yields['L'] += t13->GetEntries(PeakingCut && Cuts && TCut(LowSidebandCut.c_str()))/KpiSettings::Get().GetDouble("D0D0LuminosityScale");
  Yields['S'] += t13->GetEntries(PeakingCut && Cuts && TCut(SignalRegionCut.c_str()))/KpiSettings::Get().GetDouble("D0D0LuminosityScale");
  Yields['H'] += t13->GetEntries(PeakingCut && Cuts && TCut(UpperSidebandCut.c_str()))/KpiSettings::Get().GetDouble("D0D0LuminosityScale");
  std::cout << "Yield in region L: " << Yields['L'] << "\n";
  std::cout << "Yield in region S: " << Yields['S'] << "\n";
  std::cout << "Yield in region H: " << Yields['H'] << "\n";
  std::ifstream DataFiles(KpiSettings::Get().GetString("DataFiles"));
  TChain Chain(KpiSettings::Get().GetString("TreeName").c_str());
  while(std::getline(DataFiles, line) && KpiSettings::Get().GetString("PlotData") == "True") {
    Chain.Add(line.c_str());
  }
  DataFiles.close();
  TH1D h14("h14", "Data", 100, KpiSettings::Get().GetDouble("Min"), KpiSettings::Get().GetDouble("Max"));
  Chain.Draw((PlotVariable + std::string(" >> h14")).c_str(), Cuts, "goff");
  h1.SetLineStyle(1);
  h3.SetLineStyle(1);
  h5.SetLineStyle(1);
  h7.SetLineStyle(1);
  h9.SetLineStyle(1);
  h11.SetLineStyle(1);
  h13.SetLineStyle(1);
  h14.SetLineStyle(1);
  h1.SetLineColor(kRed);
  h3.SetLineColor(kCyan + 2);
  h5.SetLineColor(kMagenta + 4);
  h7.SetLineColor(kMagenta);
  h9.SetLineColor(kYellow - 9);
  h11.SetLineColor(kGreen);
  h13.SetLineColor(kBlue);
  h14.SetLineColor(kBlack);
  h1.SetFillColor(kRed);
  h3.SetFillColor(kCyan + 2);
  h5.SetFillColor(kMagenta + 4);
  h7.SetFillColor(kMagenta);
  h9.SetFillColor(kYellow - 9);
  h11.SetFillColor(kGreen);
  h13.SetFillColor(kBlue);
  THStack hs;
  hs.Add(&h1);
  hs.Add(&h3);
  hs.Add(&h5);
  hs.Add(&h7);
  hs.Add(&h9);
  hs.Add(&h11);
  hs.Add(&h13);
  TCanvas c1("c1", "c1", 1200, 900);
  gStyle->SetOptStat(0);
  hs.SetTitle((KpiSettings::Get().GetString("TagModeTitle") + std::string(";") + KpiSettings::Get().GetString("xLabel") + std::string(" (GeV);Events")).c_str());
  hs.SetMaximum(KpiSettings::Get().GetDouble("Height1"));
  hs.Draw("HIST");
  h14.Draw("P E1 SAME");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  c1.SaveAs(KpiSettings::Get().GetString("PlotFilename1").c_str());
  hs.SetMaximum(KpiSettings::Get().GetDouble("Height2"));
  hs.Draw("HIST");
  h14.Draw("P E1 SAME");
  gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
  c1.SaveAs(KpiSettings::Get().GetString("PlotFilename2").c_str());
  return 0;
}
