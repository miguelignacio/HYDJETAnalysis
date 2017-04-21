
import ROOT
from ROOT import TFile

import ROOT
from ROOT import TH1F
from ROOT import TCanvas
from ROOT import TLatex
import numpy as np
from AtlasCommonUtils import *


from iminuit import Minuit, describe, Struct
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.labelsize']  = 14
mpl.rcParams['lines.markersize'] = 5
from math import cos
from math import sin
from math import pi

def FromTH1toArray(histo):
  #print ' Transforming histogram ' 
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

rfile = ROOT.TFile("systematics1.5-2.0rebinX2bg.root","READ")
rfile.Print()

###2D distributions
InPlane = rfile.Get("inSignalPlusBackgroundCopy2") 
MidPlane = rfile.Get("midSignalPlusBackgroundCopy2")
OutPlane = rfile.Get("outSignalPlusBackgroundCopy2")
AllAngles = rfile.Get("allSignalPlusBackgroundCopy2")

data = {}
data["SignalRegion"] = {}
data["SignalRegion"]["AllAngles"] = FromTH1toArray(AllAngles.ProjectionX("AllAngles",24,28))
data["SignalRegion"]["InPlane"]   = FromTH1toArray(InPlane.ProjectionX("InPlane",24,28))
data["SignalRegion"]["MidPlane"]  = FromTH1toArray(MidPlane.ProjectionX("MidPlane",24,28))
data["SignalRegion"]["OutPlane"]  = FromTH1toArray(OutPlane.ProjectionX("OutPlane",24,28))


f, axes = plt.subplots(1,4, figsize=(12,3), sharex=True, sharey=True)

for counter, key in enumerate(data["SignalRegion"].keys()):
    x = data["SignalRegion"][key]["x_center"]
    y = data["SignalRegion"][key]["y"]
    dy = data["SignalRegion"][key]["dy"]
    axes[counter].set_title(key)
    axes[counter].set_xlabel(r' $\Delta\phi/\pi$')
    axes[counter].set_ylabel(r' $1/N_{trig}$ $dN/\Delta\phi$')
    axes[counter].errorbar(x,y,yerr=dy, fmt='ro')

f.subplots_adjust(wspace=0)
plt.tight_layout()
plt.show()
