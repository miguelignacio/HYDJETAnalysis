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




Correlation_Raw = rfile.Get("Correlation_Raw")
Correlation_Corrected = rfile.Get("Correlation_Corrected")


c2 = ROOT.TCanvas()
Correlation_Raw.Draw("lego2z")
Correlation_Raw.GetXaxis().SetNdivisions(4)
Correlation_Raw.GetYaxis().SetNdivisions(4)
Correlation_Raw.GetZaxis().SetNdivisions(4)
c2.SaveAs("Corr_Raw.png")
c2.Clear()
Correlation_Corrected.Draw("lego2z")
Correlation_Corrected.GetXaxis().SetNdivisions(4)
Correlation_Corrected.GetYaxis().SetNdivisions(4)
Correlation_Corrected.GetZaxis().SetNdivisions(4)
c2.SaveAs("Corr_Corr.png")




Sparse = rfile.Get("Sparse")
#defining axis
axis_dphi = 0
axis_deta = 1
axis_phitrigger = 2
axis_signal = 3
axis_pTassoc = 4

dphi_deta = {}
dphi = {}
phi_trigger = {}
deta = {}
signal = {}
pT_assoc = {}

#Sparse.GetAxis(axis_phitrigger).SetRange(1, 1) #require first bin in dphi wtr to reaction plane
#Sparse.GetAxis(axis_signal).SetRange(1,1) #require signal only
#Sparse.GetAxis(axis_deta).SetRange(11,100) #this means deta>1.
llaves_eta = {"LowEta", "LargeEta", "All"}
llaves_plane = {"Inplane", "Midplane", "Outplane", "AllAngles"}
llaves_signal = {"Signal","Bkg", "All"}

for key in llaves_signal:
    dphi_deta[key] = {}
    dphi[key] = {}
    phi_trigger[key] = {}
    deta[key] = {}
    signal[key] ={}
    pT_assoc[key] = {}

for key in llaves_signal:
    for key_eta in llaves_eta:
        dphi_deta[key][key_eta] = {}
        dphi[key][key_eta] = {}
        phi_trigger[key][key_eta] = {}
        deta[key][key_eta] = {}
        signal[key][key_eta] ={}
        pT_assoc[key][key_eta] = {}

#######################################################################################
for key_signal in llaves_signal:
    if key_signal=="Signal"       : Sparse.GetAxis(axis_signal).SetRange(2,2)
    elif key_signal=="Bkg"        : Sparse.GetAxis(axis_signal).SetRange(1,1)
    elif key_signal=="All"        : Sparse.GetAxis(axis_signal).SetRange(1,2)

    for key_plane in llaves_plane:
        if   key_plane=="Inplane"   : Sparse.GetAxis(axis_phitrigger).SetRange(1,1)
        elif key_plane=="Midplane"  : Sparse.GetAxis(axis_phitrigger).SetRange(2,2)
        elif key_plane=="Outplane"  : Sparse.GetAxis(axis_phitrigger).SetRange(3,3)
        elif key_plane=="AllAngles" : Sparse.GetAxis(axis_phitrigger).SetRange(1,3)
        for key_eta in llaves_eta:
            if   key_eta=="LowEta":   Sparse.GetAxis(axis_deta).SetRange(0,5) #this means deta < 0.5
            elif key_eta=="LargeEta": Sparse.GetAxis(axis_deta).SetRange(11,100) 
            elif key_eta=="All":      Sparse.GetAxis(axis_deta).SetRange(0,20)

            dphi_deta[key_signal][key_eta][key_plane]   = Sparse.Projection(0,1)
            phi_trigger[key_signal][key_eta][key_plane] = Sparse.Projection(2)
            dphi[key_signal][key_eta][key_plane]        = Sparse.Projection(axis_dphi)
            deta[key_signal][key_eta][key_plane]        = Sparse.Projection(axis_deta)
            signal[key_signal][key_eta][key_plane]      = Sparse.Projection(axis_signal)
            pT_assoc[key_signal][key_eta][key_plane]    = Sparse.Projection(axis_pTassoc)


def QA():

  for key_signal in llaves_signal:
    for key_plane in llaves_plane:
      for key_eta in llaves_eta:
        c1 = ROOT.TCanvas()
        c1.Divide(3,2)
        c1.cd(1)
        dphi_deta[key_signal][key_eta][key_plane].Draw("colz")
        
        c1.cd(2)
        phi_trigger[key_signal][key_eta][key_plane].SetTitle("; phi_trigger; entries")
        phi_trigger[key_signal][key_eta][key_plane].Draw()
        phi_trigger[key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(6)
        
        c1.cd(3)
        dphi[key_signal][key_eta][key_plane].SetLineColor(2)
        dphi[key_signal][key_eta][key_plane].SetTitle("; #Delta#phi/#pi; entries")
        dphi[key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(6)
        dphi[key_signal][key_eta][key_plane].GetYaxis().SetNdivisions(4)
        dphi[key_signal][key_eta][key_plane].Draw()
        c1.cd(4)
        deta[key_signal][key_eta][key_plane].SetTitle("; deta; entries")
        deta[key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(4)
        deta[key_signal][key_eta][key_plane].Draw()
        
        c1.cd(5)
        signal[key_signal][key_eta][key_plane].SetTitle("; signal or not ; entries")
        signal[key_signal][key_eta][key_plane].Draw()
        signal[key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(2)
        c1.cd(6)
        pT_assoc[key_signal][key_eta][key_plane].SetTitle("; pT assoc. [GeV] ; entries")
        pT_assoc[key_signal][key_eta][key_plane].Draw()
        pT_assoc[key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(5)
        c1.SaveAs("QA_%s_%s_%s.png" %(key_signal,key_eta,key_plane))
        return


def DrawHistos():
  for key_eta in llaves_eta:
    c = ROOT.TCanvas()
    c.Divide(4,3)
    c.cd(1)
    dphi["All"][key_eta]["Inplane"].Draw()
    dphi["All"][key_eta]["Inplane"].SetLineColor(1)
    c.cd(2)
    dphi["All"][key_eta]["Midplane"].Draw()
    dphi["All"][key_eta]["Midplane"].SetLineColor(1)
    c.cd(3)
    dphi["All"][key_eta]["Outplane"].Draw()
    dphi["All"][key_eta]["Outplane"].SetLineColor(1)
    c.cd(4)
    dphi["All"][key_eta]["AllAngles"].Draw()
    dphi["All"][key_eta]["AllAngles"].SetLineColor(1)
    
    c.cd(5)
    dphi["Bkg"][key_eta]["Inplane"].Draw()
    dphi["Bkg"][key_eta]["Inplane"].SetLineColor(2)
    c.cd(6)
    dphi["Bkg"][key_eta]["Midplane"].Draw()
    dphi["Bkg"][key_eta]["Midplane"].SetLineColor(2)
    c.cd(7)
    dphi["Bkg"][key_eta]["Outplane"].Draw()
    dphi["Bkg"][key_eta]["Outplane"].SetLineColor(2)
    c.cd(8)
    dphi["Bkg"][key_eta]["AllAngles"].Draw()
    dphi["Bkg"][key_eta]["AllAngles"].SetLineColor(2)
    
    c.cd(9)
    dphi["Signal"][key_eta]["Inplane"].Draw()
    dphi["Signal"][key_eta]["Inplane"].SetLineColor(4)
    c.cd(10)
    dphi["Signal"][key_eta]["Midplane"].Draw()
    dphi["Signal"][key_eta]["Midplane"].SetLineColor(4)
    c.cd(11)
    dphi["Signal"][key_eta]["Outplane"].Draw()
    dphi["Signal"][key_eta]["Outplane"].SetLineColor(4)
    c.cd(12)
    dphi["Signal"][key_eta]["AllAngles"].Draw()
    dphi["Signal"][key_eta]["AllAngles"].SetLineColor(4)
    c.SaveAs("SUMMARY_%s.png" %(key_eta))
            

def FromTH1toArray(histo):
  print ' Transforming histogram ' 
  # Calculate how many y bins we have
  xaxis = histo.GetXaxis()
  nx = xaxis.GetNbins()
  xbins = range(1, nx+1) # Removes under/overflow
  d = {}
  d["x_center"] = [xaxis.GetBinCenter(i) for i in xbins]
  d["x_low"]    = [xaxis.GetBinLowEdge(i) for i in xbins]
  d["x_up"]     = [xaxis.GetBinUpEdge(i) for i in xbins]
  d["y"]        = [histo.GetBinContent(i) for i in xbins]
  d['dy']       = [histo.GetBinError(i) for i in xbins]
  return d
   

QA()
#datos = FromTH1toArray(dphi["Bkg"][key_eta]["Midplane"])
#import matplotlib.pyplot as plt
#plt.plot(datos["x_center"], datos["y"], '-ro')
#plt.show()
