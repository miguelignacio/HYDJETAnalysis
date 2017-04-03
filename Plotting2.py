import ROOT
from ROOT import TH1F
from ROOT import TCanvas
from ROOT import TLatex
import numpy as np
from AtlasCommonUtils import *


from iminuit import Minuit
import matplotlib.pyplot as plt
from math import cos



def func(x, a,b):
    return a*(1 + b*np.cos(2*np.pi*x))


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
  d["x_center"] = [xaxis.GetBinCenter(i) for i in xbins]
  d["x_low"]    = [xaxis.GetBinLowEdge(i) for i in xbins]
  d["x_up"]     = [xaxis.GetBinUpEdge(i) for i in xbins]
  d["y"]        = [histo.GetBinContent(i) for i in xbins]
  d['dy']       = [histo.GetBinError(i) for i in xbins]
  return d
   

#QA()
#DrawHistos(histo["dphi"])

#datos = {}
#datos["All"] = {}
#datos["Bkg"] = {}
#datos["All"]["LowEta"] = FromTH1toArray(histo["dphi"]["All"]["LowEta"]["AllAngles"])
#datos["Bkg"]["LowEta"] = FromTH1toArray(histo["dphi"]["Bkg"]["LowEta"]["AllAngles"])
#datos["All"]["LargeEta"] = FromTH1toArray(histo["dphi"]["All"]["LargeEta"]["AllAngles"])
#datos["Bkg"]["LargeEta"] = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["AllAngles"])



# Two subplots, unpack the axes array immediately

from math import cos

print ' COS ' , cos(0)

#def FitFunction(x, v2_trig, v3_trig, v2_assoc, v3_assoc, B):
#def FitFunction(x,v2_trig, B):
#  return 7
  #return B*(1+ v2_trig*cos(2*x))
    #temp = B* (1+ v2_trig*v2_assoc*cos(2*x) + v3_trig*v3_assoc*cos(3*x))
 #   return temp


#def ModV(n, phi_s, c, R):
#  primer  = v[n] + R*cos(n*phi_s)*sin(n*c)/(nc)
#  segundo = R*(v[abs(k+n)]+v[abs(k-n)])*cos(k*phi_s)*sin(k*c)/(k*c)
#  tercer  = 1 + 2*R*v[k]*cos(k*phi_s)*sin(k*c)/(k*c)


from scipy.optimize import curve_fit


h1 = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["AllAngles"])
h2 = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["Inplane"])
h3 = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["Midplane"])
h4 = FromTH1toArray(histo["dphi"]["Bkg"]["LargeEta"]["Outplane"])

x = np.array(h2["x_center"])
y = np.array(h2["y"])

print x
print y


def MyCHI(a,b):
    return np.sum( np.power(y -func(x, a,b)  , 2.0))

#m = Minuit(MyCHI, a=0.8, b=0.8, limit_a=(-1,3), limit_b=(-1,5))

m = Minuit(MyCHI)
m.migrad()
print ' values'
print(m.values)  # {'x': 2,'y': 3,'z': 4}
print ' errors'
print(m.errors)  # {'x': 1,'y': 1,'z': 1}



plt.plot(x, func(x, m.values["a"], m.values["b"]) , "-")
#plt.plot(x, y, 'o')
#plt.show()


f, (ax1, ax2, ax3) = plt.subplots(1,3, sharex=True, sharey=True)
ax1.plot(h2["x_center"], h2["y"], '-o')
ax1.plot(x, func(x, m.values["a"], m.values["b"]) , "-")
ax2.plot(h3["x_center"], h3["y"], '-o')
ax3.plot(h4["x_center"], h4["y"], '-o')





#for key_planes in llaves_plane:
#  h = FromTH1toArray(histo["dphi"]["All"]["LowEta"][key_planes])
#  plt.plot(h["x_center"], h["y"], "-o")


#plt.plot(datos["All"]["LowEta"]["x_center"], datos["All"]["LowEta"]["y"], '-ro')
#plt.plot(datos["All"]["LargeEta"]["x_center"], datos["All"]["LargeEta"]["y"], '-ro')
plt.show()
