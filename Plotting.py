import ROOT
from ROOT import TH1F
from ROOT import TCanvas
from ROOT import TLatex
import numpy as np
from AtlasCommonUtils import *

def SubtractPedestal(h):
  minimo = h.GetMinimum();
  print minimo
  for n in range(h.GetNbinsX()+1):
      h.SetBinContent(n, (h.GetBinContent(n)-minimo)/1.0);
      
  h.Scale(1.0/1000.0)  
  return

SetAtlasStyle()

rfile = ROOT.TFile("fout_histos.root","READ")
rfile.Print()

h_dphi_signal = rfile.Get("Dphi_Signal")
h_dphi_all =    rfile.Get("Dphi_All")
h_dphi_bkg = rfile.Get("Dphi_Background")



c1 = ROOT.TCanvas("c1","c1",600,600)
#SubtractPedestal(h_dphi_signal)
#h_dphi_signal.Draw()
SubtractPedestal(h_dphi_all)
h_dphi_all.Draw()
h_dphi_all.SetLineColor(4)
h_dphi_all.SetMarkerColor(4)
h_dphi_all.GetXaxis().SetNdivisions(5)
h_dphi_all.GetYaxis().SetNdivisions(5)
h_dphi_all.SetTitle("; #Delta#phi/#pi; C(#Delta#phi)")
SubtractPedestal(h_dphi_bkg)
h_dphi_bkg.SetLineColor(2)
h_dphi_bkg.SetMarkerColor(2)

h_dphi_all.SetMarkerStyle(20)
h_dphi_bkg.SetMarkerStyle(20)
h_dphi_all.SetMarkerSize(1)
h_dphi_bkg.SetMarkerSize(1)

latex = TLatex();
latex.SetNDC();
latex.DrawLatex(0.52,0.56, "#font[42]{#scale[0.9]{#color[2]{Background }}}");
latex.DrawLatex(0.52,0.62, "#font[42]{#scale[0.9]{#color[4]{Signal + Background }}}");
latex.DrawLatex(0.15,0.96, "#font[42]{#scale[0.8]{30-40% Pb+Pb 2.7 TeV,}}")
latex.DrawLatex(0.55,0.96, "#font[42]{#scale[0.8]{PYTHIA + HYDJET}}")
latex.DrawLatex(0.52,0.77, "#font[42]{#scale[0.9]{%2.f < p_{T}^{b} < %2.f GeV}}" %(1.0,2.0))
latex.DrawLatex(0.52,0.84, "#font[42]{#scale[0.9]{%2.f < p_{T}^{a} < %2.f GeV}}" %(6.0,8.0\
))
latex.DrawLatex(0.57,0.71, "#font[42]{#scale[0.9]{|#Delta#eta| < %3.2f}}" %(0.5))
h_dphi_bkg.Draw("same")






c1.SaveAs("histogram.png")
