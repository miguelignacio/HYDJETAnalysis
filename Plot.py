import ROOT
from ROOT import TH1F
from ROOT import TCanvas
from ROOT import TLatex
import numpy as np
from AtlasCommonUtils import *


from iminuit import Minuit
import matplotlib.pyplot as plt
from math import cos



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


Sparse = rfile.Get("Sparse")
Sparse_trigger = rfile.Get("Sparse_trigger")

#defining axis
axis_dphi = 0
axis_deta = 1
axis_phitrigger = 2
axis_signal = 3
axis_pTassoc = 4
axis_pTtrig  = 5


histo = {}

histo["dphi_deta"]   = {}
histo["dphi"]        = {}
histo["phi_trigger"] = {}
histo["deta"]        = {}
histo["signal"]      = {}
histo["pT_assoc"]    = {}
histo["pT_trig"]     = {}


#Sparse.GetAxis(axis_phitrigger).SetRange(1, 1) #require first bin in dphi wtr to reaction plane
#Sparse.GetAxis(axis_signal).SetRange(1,1) #require signal only
#Sparse.GetAxis(axis_deta).SetRange(11,100) #this means deta>1.
Sparse.GetAxis(axis_pTtrig).SetRange(4,100) #pT trigger > 5 GeV
#Sparse.GetAxis(0).SetRange(6,100)

Sparse.GetAxis(axis_pTassoc).SetRange(3,4) # pT associated 1--2 GeV
#llaves_eta = {"LowEta", "LargeEta", "All"}
llaves_eta = {"LargeEta"}
llaves_plane = {"Inplane", "Midplane", "Outplane", "AllAngles"}
#llaves_signal = {"Signal","Bkg", "All"}
llaves_signal = {"Bkg"}


nTriggers = {}
nTriggers["AllAngles"] = Sparse_trigger.Projection(0).GetEntries()
nTriggers["Inplane"]   = Sparse_trigger.Projection(1).GetBinContent(1)
nTriggers["Midplane"]  = Sparse_trigger.Projection(1).GetBinContent(2)
nTriggers["Outplane"]  = Sparse_trigger.Projection(1).GetBinContent(3)

for key in nTriggers.keys():
  print key, nTriggers[key]


for key in llaves_signal:
    histo["dphi_deta"][key]   = {}
    histo["dphi"][key]        = {}
    histo["phi_trigger"][key] = {}
    histo["deta"][key]        = {}
    histo["signal"][key]      = {}
    histo["pT_assoc"][key]    = {}
    histo["pT_trig"][key]     = {}

for key in llaves_signal:
    for key_eta in llaves_eta:
        histo["dphi_deta"][key][key_eta] = {}
        histo["dphi"][key][key_eta] = {}
        histo["phi_trigger"][key][key_eta] = {}
        histo["deta"][key][key_eta] = {}
        histo["signal"][key][key_eta] ={}
        histo["pT_assoc"][key][key_eta] = {}
        histo["pT_trig"][key][key_eta] = {}

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

            histo["dphi_deta"][key_signal][key_eta][key_plane]   = Sparse.Projection(axis_dphi, axis_deta)
            histo["phi_trigger"][key_signal][key_eta][key_plane] = Sparse.Projection(axis_phitrigger)
            histo["dphi"][key_signal][key_eta][key_plane]        = Sparse.Projection(axis_dphi)
            histo["dphi"][key_signal][key_eta][key_plane].Scale(1.0/nTriggers[key_plane])

            histo["deta"][key_signal][key_eta][key_plane]        = Sparse.Projection(axis_deta)
            histo["signal"][key_signal][key_eta][key_plane]      = Sparse.Projection(axis_signal)
            histo["pT_assoc"][key_signal][key_eta][key_plane]    = Sparse.Projection(axis_pTassoc)
            histo["pT_trig"][key_signal][key_eta][key_plane]    =  Sparse.Projection(axis_pTtrig)
            


def QA():

  llaves_eta = {"LowEta", "LargeEta", "All"}
  llaves_plane = {"Inplane", "Midplane", "Outplane", "AllAngles"}
  llaves_signal = {"Signal","Bkg", "All"}
  
  for key_signal in llaves_signal:
    for key_plane in llaves_plane:
      for key_eta in llaves_eta:
        print " %s %s %s " %(key_signal,key_plane,key_eta)
        c1 = ROOT.TCanvas()
        c1.Divide(3,3)
        c1.cd(1)
        histo["dphi_deta"][key_signal][key_eta][key_plane].Draw("colz")
        
        c1.cd(2)
        histo["phi_trigger"][key_signal][key_eta][key_plane].SetTitle("; phi_trigger; entries")
        histo["phi_trigger"][key_signal][key_eta][key_plane].Draw()
        histo["phi_trigger"][key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(6)
        
        c1.cd(3)
        histo["dphi"][key_signal][key_eta][key_plane].SetLineColor(2)
        histo["dphi"][key_signal][key_eta][key_plane].SetTitle("; #Delta#phi/#pi; entries")
        histo["dphi"][key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(6)
        histo["dphi"][key_signal][key_eta][key_plane].GetYaxis().SetNdivisions(4)
        histo["dphi"][key_signal][key_eta][key_plane].Draw()
        c1.cd(4)
        histo["deta"][key_signal][key_eta][key_plane].SetTitle("; deta; entries")
        histo["deta"][key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(4)
        histo["deta"][key_signal][key_eta][key_plane].Draw()
        
        c1.cd(5)
        histo["signal"][key_signal][key_eta][key_plane].SetTitle("; signal or not ; entries")
        histo["signal"][key_signal][key_eta][key_plane].Draw()
        histo["signal"][key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(2)
        c1.cd(6)
        histo["pT_assoc"][key_signal][key_eta][key_plane].SetTitle("; pT assoc. [GeV] ; entries")
        histo["pT_assoc"][key_signal][key_eta][key_plane].Draw()
        histo["pT_assoc"][key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(5)

        c1.cd(7)
        histo["pT_trig"][key_signal][key_eta][key_plane].SetTitle("; pT trig [GeV] ; entries")
        histo["pT_trig"][key_signal][key_eta][key_plane].Draw()
        histo["pT_trig"][key_signal][key_eta][key_plane].GetXaxis().SetNdivisions(5)
        c1.SaveAs("QA_%s_%s_%s.pdf" %(key_signal,key_eta,key_plane))



colores = {}
colores["All"] = 1
colores["Signal"] = 4
colores["Bkg"] = 2 
def DrawHistos(histo,option=""):
  for key_eta in llaves_eta:
    c = ROOT.TCanvas()
    c.Divide(4,3)
    it = 1
    #for key_signal in llaves_signal:
    #key_signal = "Signal"
    for key_signal in llaves_signal:
      for key_plane in llaves_plane:
        c.cd(it)
        h = histo[key_signal][key_eta][key_plane]
        h.Draw("hist")
        h.SetLineColor(colores[key_signal])
        h.SetMarkerColor(colores[key_signal])
        h.Scale(1.0/nTriggers[key_plane])
        it = it+1
    c.SaveAs("SUMMARY_%s.pdf" %(key_eta))

def DrawHistos2(histo,option=""):
  for key_eta in llaves_eta:
    c = ROOT.TCanvas()
    c.Divide(4,2)
    it = 1    
    #for key_plane in llaves_plane:
    #  c.cd(it)
    hs = ROOT.THStack("hs","hs")
    c.cd(it)
    key_plane = "AllAngles"
    histo["All"][key_eta][key_plane].Scale(1.0/nTriggers[key_plane])
    histo["Bkg"][key_eta][key_plane].Scale(1.0/nTriggers[key_plane])
    histo["Bkg"][key_eta][key_plane].SetFillColor(5)
    histo["Bkg"][key_eta][key_plane].SetLineColor(4)
    hs.Add(histo["Signal"][key_eta][key_plane])
    hs.Add(histo["Bkg"][key_eta][key_plane])
    
    #h = histo[key_signal][key_eta][key_plane]
    histo["All"][key_eta][key_plane].Draw()
    histo["Bkg"][key_eta][key_plane].Draw("same")
    #hs.Draw("nostack")
    #h.SetLineColor(colores[key_signal])
    #h.SetMarkerColor(colores[key_signal])
    #h.Scale(1.0/nTriggers[key_plane])
    it = it+1
    c.SaveAs("SUMMARY_%s.pdf" %(key_eta))
    

def FromTH1toArray(histo):
  print ' Transforming histogram ' 
  # Calculate how many y bins we have
  xaxis = histo.GetXaxis()
  nx = xaxis.GetNbins()
  xbins = range(1, nx+1) # Removes under/overflow
  d = {}
  d["x_center"] = np.array([xaxis.GetBinCenter(i) for i in xbins])
  d["x_low"]    = np.array([xaxis.GetBinLowEdge(i) for i in xbins])
  d["x_up"]     = np.array([xaxis.GetBinUpEdge(i) for i in xbins])
  d["y"]        = np.array([histo.GetBinContent(i) for i in xbins])
  d['dy']       = np.array([histo.GetBinError(i) for i in xbins])
  return d
   

#QA()
#DrawHistos(histo["dphi"])

# Two subplots, unpack the axes array immediately

from math import cos

print ' COS ' , cos(0)

from scipy.optimize import curve_fit


data = {}

#data["AllAngles"] = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["AllAngles"])
data["Inplane"]   = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["Inplane"])
data["Midplane"]  = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["Midplane"])
data["Outplane"]  = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["Outplane"])

def func(x, B, v2_trig, v2_assoc, V3, v4_trig, v4_assoc, phi_s, c):

    num = v2_trig + np.cos(2*phi_s)*np.sin(2*c)/(2*c) +  v4_trig*np.cos(2*phi_s)*np.sin(2*c)/(2*c) +  v2_trig*np.cos(4*phi_s)*np.sin(4*c)/(4*c) +   v4_trig*np.cos(6*phi_s)*np.sin(6*c)/(6*c)

    den =  1 + 2*v2_trig*np.cos(2*phi_s)*np.sin(2*c)/(2*c) + 2*v4_trig*np.cos(4*phi_s)*np.sin(4*c)/(2*c) 
    v2R = num/den
	
    num2 = v4_trig +  np.cos(4*phi_s)*np.sin(4*c)/(4*c) +    v2_trig*np.cos(2*phi_s)*np.sin(2*c)/(2*c) +   v2_trig*np.cos(6*phi_s)*np.sin(6*c)/(6*c) +  v4_trig*np.cos(8*phi_s)*np.sin(8*c)/(8*c) 
	
    v4R = num2/den
	
    #B = B*den

    return B*(1 + v2R*v2_assoc*np.cos(2*np.pi*x) + V3*np.cos(3*np.pi*x) + v4R*v4_assoc*np.cos(4*np.pi*x) )

def MyCHI(B, v2_trig, v2_assoc, V3, v4_trig, v4_assoc):

    total_chi2 = 0
    phi_s = {}
    phi_s["Inplane"]  = 0
    phi_s["Midplane"] =np.pi/4.0
    phi_s["Outplane"] = np.pi/2.0
    c  = {}
    c["Inplane"]  = np.pi/6.0 
    c["Midplane"] = np.pi/12.0
    c["Outplane"] = np.pi/6.0
    for key in data.keys():
        if key=="Midplane": continue
        x = data["Inplane"]["x_center"]
        y = data["Inplane"]["y"]
        
        model = func( x , B, v2_trig, v2_assoc, V3, v4_trig, v4_assoc, phi_s[key], c[key] )
       
        total_chi2 = total_chi2+ np.sum( np.power( y - model  , 2.0))
    return total_chi2

#m = Minuit(MyCHI, a=0.8, b=0.8, limit_a=(-1,3), limit_b=(-1,5))


print 'About to start MINUIT ' 

m = Minuit(MyCHI)
m.migrad()
print ' values'
print(m.values)  # {'x': 2,'y': 3,'z': 4}
print ' errors'
print(m.errors)  # {'x': 1,'y': 1,'z': 1}



f, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(1,6, sharex=True, sharey=True)
ax1.plot(data["Inplane"]["x_center"],  data["Inplane"]["y"], '-o')
ax2.plot(data["Midplane"]["x_center"], data["Midplane"]["y"], '-o')
ax3.plot(data["Outplane"]["x_center"], data["Outplane"]["y"], '-o')

pre1 = func(data["Inplane"]["x_center"], m.values["B"], m.values["v2_trig"], m.values["v2_assoc"], m.values["V3"], 
	m.values["v4_trig"], m.values["v4_assoc"], 0, np.pi/6.0) 
pre2 = func(data["Midplane"]["x_center"], m.values["B"], m.values["v2_trig"], m.values["v2_assoc"], m.values["V3"], 
	m.values["v4_trig"], m.values["v4_assoc"], 	np.pi/4.0, np.pi/12.0)
pre3 = func(data["Outplane"]["x_center"], m.values["B"], m.values["v2_trig"], m.values["v2_assoc"], m.values["V3"],
	m.values["v4_trig"], m.values["v4_assoc"], np.pi/2.0, np.pi/6.0)

ax4.plot(data["Inplane"]["x_center"],  pre1, '-')
ax5.plot(data["Midplane"]["x_center"],  pre2, '-')
ax6.plot(data["Outplane"]["x_center"],  pre3, '-')


#for key_planes in llaves_plane:
#  h = FromTH1toArray(histo["dphi"]["All"]["LowEta"][key_planes])
#  plt.plot(h["x_center"], h["y"], "-o")


#plt.plot(datos["All"]["LowEta"]["x_center"], datos["All"]["LowEta"]["y"], '-ro')
#plt.plot(datos["All"]["LargeEta"]["x_center"], datos["All"]["LargeEta"]["y"], '-ro')
plt.show()
