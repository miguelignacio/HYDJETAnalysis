import ROOT
from ROOT import TH1F
from ROOT import TCanvas
from ROOT import TLatex
import numpy as np
from AtlasCommonUtils import *


from iminuit import Minuit, describe, Struct
import matplotlib.pyplot as plt

import matplotlib as mpl
from matplotlib import colors as mcolors

mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['axes.labelsize']  = 14
#mpl.rcParams['legend.fontsize'] = 0.5 
#mpl.rcParams['legend.linewidth'] = 0.8
mpl.rcParams['lines.markersize'] = 3.8
mpl.rcParams['legend.numpoints'] =1
mpl.rcParams['legend.fontsize'] = 'small'
from math import cos
from math import sin
from math import pi


from ROOT import gROOT
gROOT.ProcessLine( "gErrorIgnoreLevel = 2001;") #asdasd
print 'Ola Mundo'



phi_s = {}
phi_s["Inplane"]  = 0
phi_s["Midplane"] =np.pi/4.0
phi_s["Outplane"] = np.pi/2.0
#phi_s["AllAngles"]  = 0
c  = {}
c["Inplane"]  = np.pi/6.0
c["Midplane"] = np.pi/12.0
c["Outplane"] = np.pi/6.0
#c["AllAngles"] = 2*np.pi


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


def SubtractPedestal(h):
  minimo = h.GetMinimum();
  print minimo
  for n in range(h.GetNbinsX()+1):
      h.SetBinContent(n, (h.GetBinContent(n)-minimo)/1.0);
      
  h.Scale(1.0/1000.0)  
  return

SetAtlasStyle()

rfile  = ROOT.TFile("ROOTFILES/fout_Cen3040_20k_Baseline.root","READ")
rfile.Print()
Sparse = {}
Sparse_trigger = {}
Sparse["Baseline"]  = rfile.Get("Sparse")
Sparse_trigger["Baseline"] = rfile.Get("Sparse_trigger")


rfile = ROOT.TFile("ROOTFILES/fout_Cen3040_20k_Baseline.root","READ")
rfile.Print()
Sparse["Quenched"]  = rfile.Get("Sparse")
Sparse_trigger["Quenched"] = rfile.Get("Sparse_trigger")

#defining axis
axis_dphi = 0
axis_deta = 1
axis_phitrigger = 2
axis_signal = 3
axis_pTassoc = 4
axis_pTtrig  = 5


#llaves_eta = {"LowEta", "LargeEta", "All"}
llaves_eta = {"LowEta", "LargeEta"}
#llaves_eta = {"LargeEta"}
#llaves_plane = {"Inplane", "Midplane", "Outplane", "AllAngles"}
llaves_plane = {"Inplane", "Midplane", "Outplane"}
llaves_signal = {"Signal","Bkg", "All"}
#llaves_signal = {"Bkg"}

def GettingDataFromFile(Sparse, Sparse_trigger):

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
    ptmin = 6
    ptmax = 10
    Sparse.GetAxis(axis_pTtrig).SetRange(ptmin+1,ptmax) #pT trigger  5--10 GeV
    Sparse_trigger.GetAxis(0).SetRange(ptmin+1, ptmax) ##the same for the histogram

    Sparse.GetAxis(axis_pTassoc).SetRange(3,4) # pT associated 1--2 GeV (3, 4)

    nTriggers = {}
    nTriggers["AllAngles"] = Sparse_trigger.Projection(0).GetEntries()
    nTriggers["Inplane"]   = Sparse_trigger.Projection(1).GetBinContent(1)
    nTriggers["Midplane"]  = Sparse_trigger.Projection(1).GetBinContent(2)
    nTriggers["Outplane"]  = Sparse_trigger.Projection(1).GetBinContent(3)

    for key in nTriggers.keys():
        print 'Number of Triggers' , key, ' ' , nTriggers[key]

    #for key_signal in llaves_signal:
    #  for htype in histo.keys():
    #      histo[htype][key_signal]   = {}
    #Initializing histograms:
    for key_signal in llaves_signal:
        for key_type in histo.keys():
            histo[key_type][key_signal] = {}
            for key_eta in llaves_eta:
                histo[key_type][key_signal][key_eta] = {}

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
                if   key_eta=="LowEta":   Sparse.GetAxis(axis_deta).SetRange(0,5) # (0,5) this means deta < 0.5
                elif key_eta=="LargeEta": Sparse.GetAxis(axis_deta).SetRange(11,15) #(11,15) this means from 1.0 --1.5  
                elif key_eta=="All":      Sparse.GetAxis(axis_deta).SetRange(0,20)

                histo["dphi_deta"][key_signal][key_eta][key_plane]   = Sparse.Projection(axis_dphi, axis_deta)
                histo["phi_trigger"][key_signal][key_eta][key_plane] = Sparse.Projection(axis_phitrigger)
                histo["dphi"][key_signal][key_eta][key_plane]        = Sparse.Projection(axis_dphi)

                histo["deta"][key_signal][key_eta][key_plane]        = Sparse.Projection(axis_deta)
                histo["signal"][key_signal][key_eta][key_plane]      = Sparse.Projection(axis_signal)
                histo["pT_assoc"][key_signal][key_eta][key_plane]    = Sparse.Projection(axis_pTassoc)
                histo["pT_trig"][key_signal][key_eta][key_plane]    =  Sparse.Projection(axis_pTtrig)

                histo["dphi"][key_signal][key_eta][key_plane].Sumw2()
                #histo["dphi"][key_signal][key_eta][key_plane].Scale(1.0/nTriggers[key_plane])

    return histo        


def QA(histo):

  
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
        c1.SaveAs("PDFOUTPUT/QA_%s_%s_%s.pdf" %(key_signal,key_eta,key_plane))

    return #end of QA function





histo = {}
histo["Baseline"] = GettingDataFromFile(Sparse["Baseline"], Sparse_trigger["Baseline"])
histo["Quenched"] = GettingDataFromFile(Sparse["Quenched"], Sparse_trigger["Quenched"])


        
#QA(histo["Quenched"])

# Two subplots, unpack the axes array immediately

from math import cos

print ' COS ' , cos(0), pi

from scipy.optimize import curve_fit
from matplotlib import rc


print ' Going to start plotting part ' 


#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


def Plot(h, axes):
    for i, key_plane in enumerate(llaves_plane):
        temp = FromTH1toArray(h[key_plane])
        axes[i].errorbar(temp["x_center"],temp["y"], yerr=temp["dy"], fmt='-o')
        axes[i].locator_params(nbins=4)
        axes[i].set_title(key_plane)
        axes[i].set_xlabel(r' $\Delta\phi/\pi$')
        axes[i].set_ylabel(r' $1/N_{trig}$ $dN/\Delta\phi$')

#Here I want to plot the large-eta region for background and signal + background only.
def PlotHisto(histo):
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,3), sharex=True, sharey=True)
    Plot(histo["Bkg"]["LargeEta"], axs)
    Plot(histo["All"]["LargeEta"], axs) 
    Plot(histo["Signal"]["LargeEta"], axs)
    fig.subplots_adjust(wspace=0)
    plt.tight_layout()
    #plt.show()
    #fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,3), sharex=True, sharey=True)
    #Plot(histo["Signal"]["LargeEta"], axs)
    #fig.subplots_adjust(wspace=0)
    #plt.tight_layout()
    fig.savefig("PLOTS/LargeEta.png")
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,3), sharex=True, sharey=True)
    Plot(histo["Bkg"]["LowEta"], axs)
    Plot(histo["All"]["LowEta"], axs)
    Plot(histo["Signal"]["LowEta"], axs)
    fig.subplots_adjust(wspace=0)
    plt.tight_layout()
    #plt.show()
    fig.savefig("PLOTS/LowEta.png")
    #fig.clf()

    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,3), sharex=True, sharey=True)
    Plot(histo["Signal"]["LowEta"], axs)
    fig.subplots_adjust(wspace=0)
    plt.tight_layout()
    #plt.show()
    fig.savefig("Signal.pdf")
    #fig.clf()
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,3), sharex=True, sharey=True)
    Plot(histo["Bkg"]["LowEta"], axs)
    Plot(histo["Bkg"]["LargeEta"], axs)
    Plot(histo["All"]["LowEta"], axs)
    Plot(histo["All"]["LargeEta"], axs)
    fig.subplots_adjust(wspace=0)
    plt.tight_layout()
    #plt.show()
    fig.savefig("PLOTS/BKG.png")
    
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,3), sharex=True, sharey=True)
    Plot(histo["All"]["LargeEta"], axs)
    Plot(histo["All"]["LowEta"], axs)
    Plot(histo["Signal"]["LowEta"], axs)
    fig.subplots_adjust(wspace=0)
    plt.tight_layout()
    #plt.show()
    fig.savefig("PLOTS/SignalAndBackground.png")




def PlotHistoComparison(h1,h2):
    fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,3), sharex=True, sharey=True)
    Plot(h1["Signal"]["LowEta"], axs)
    Plot(h2["Signal"]["LowEta"], axs)
    fig.subplots_adjust(wspace=0)
    plt.tight_layout()
    plt.show()
    fig.savefig("Comparison.pdf")

#PlotHisto(histo["Quenched"]["dphi"])
#PlotHisto(histo["Baseline"]["dphi"])

#PlotHistoComparison(histo["Quenched"]["dphi"], histo["Baseline"]["dphi"])

def Comparison(histo):
    for j, key_signal in enumerate(llaves_signal):
        print j, key_signal
        fig, axs = plt.subplots(nrows=1, ncols=4, figsize=(12,3), sharex=True, sharey=True)
        for i, key_plane in enumerate(llaves_plane):
            h1 =  FromTH1toArray(histo["Baseline"]["dphi"][key_signal]["LowEta"][key_plane])
            h2 =  FromTH1toArray(histo["Quenched"]["dphi"][key_signal]["LowEta"][key_plane])
            axs[i].errorbar(h1["x_center"],h1["y"], yerr=h1["dy"], fmt='-o', label='No E-loss')
            axs[i].errorbar(h2["x_center"],h2["y"], yerr=h2["dy"], fmt='-|', label='With E-loss')
            axs[i].locator_params(nbins=4)
            axs[i].set_title(key_plane)
            axs[i].set_xlabel(r' $\Delta\phi/\pi$')
            axs[i].set_ylabel(r' $1/N_{trig}$ $dN/\Delta\phi$')
            #axs[3].legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
        fig.subplots_adjust(wspace=0)
        plt.tight_layout()
        plt.show()
        fig.savefig("PLOTS/Comparison%s.png" %(key_signal))
        #fig.savefig("PLOTS/Comparison%s.pdf" %(key_signal))
        #plt.clf()





#######################################################################################
############################### FITTING PART ##########################################
######################################################################################
dataALL = {}
dataSIGNAL = {}
dataBKG = {}
dataBKG_True = {}
#print histo["Baseline"]["dphi"]


##DEFINING DATA FOR THE FITS:
for key in llaves_plane:
    dataALL[key] = FromTH1toArray(histo["Baseline"]["dphi"]["All"]["LowEta"][key])
    dataBKG[key] = FromTH1toArray(histo["Baseline"]["dphi"]["All"]["LargeEta"][key])

    dataSIGNAL[key] = FromTH1toArray(histo["Baseline"]["dphi"]["Signal"]["LowEta"][key])
    dataBKG_True[key] =  FromTH1toArray(histo["Baseline"]["dphi"]["Bkg"]["LowEta"][key])
    print '\n'


def TotalPDF(x, params_Signal, params_BKG , phi,c):
    return Signal(x,params_Signal) + Background(x,params_BKG, phi, c)

def Signal(x, params):

    A1 = params["A1"]
    A2 = params["A2"]
    s1 = params["s1"]
    s2 = params["s2"]
    C1 = params["C1"]
    return A1*np.exp(-(x-0.0)**2/(2*s1**2) ) + A2*np.exp(-(x-1.0)**2/(2*s2**2) ) + C1

def Background(x, params, phi , c ):
    B    = params["B"]
    v2_t = params["v2_t"]
    v2_a = params["v2_a"]
    V3   = params["V3"]
    v4_t = params["v4_t"]
    v4_a = params["v4_a"]
    num = v2_t + cos(2*phi)*sin(2*c)/(2*c) +  v4_t*cos(2*phi)*sin(2*c)/(2*c) +  v2_t*cos(4*phi)*sin(4*c)/(4*c) + v4_t*cos(6*phi)*sin(6*c)/(6*c)
    den =  1 + 2*v2_t*cos(2*phi)*sin(2*c)/(2*c) + 2*v4_t*cos(4*phi)*sin(4*c)/(4*c)
    v2R = num/den
    num2 = v4_t + cos(4*phi)*sin(4*c)/(4*c) + v2_t*cos(2*phi)*sin(2*c)/(2*c)+ v2_t*cos(6*phi)*sin(6*c)/(6*c) + v4_t*cos(8*phi)*sin(8*c)/(8*c)
    v4R = num2/den
    BR = B*den*c*2/np.pi
    factor = 1.0
    if(c==np.pi/12.0): factor=2.0
    BR = BR*factor
    return BR*(1 + 2*v2R*v2_a*np.cos(2*np.pi*x) + 2*V3*np.cos(3*np.pi*x) + 2*v4R*v4_a*np.cos(4*np.pi*x))

def Chi2(A3,A4,s3,s4,C2,
         A5,A6,s5,s6,C3,
         A7,A8,s7,s8,C4,
         B, v2_t, v2_a, V3, v4_t, v4_a):
    params_Signal = {}
    params_BKG = {}

    params_Signal["Midplane"]  = {}
    params_Signal["Inplane"]  = {}
    params_Signal["Outplane"] = {}
 
    params_Signal["Outplane"]["A1"] = A3
    params_Signal["Outplane"]["A2"] = A4
    params_Signal["Outplane"]["s1"] = s3
    params_Signal["Outplane"]["s2"] = s4
    params_Signal["Outplane"]["C1"] = C2
    
    params_Signal["Midplane"]["A1"] = A5
    params_Signal["Midplane"]["A2"] = A6
    params_Signal["Midplane"]["s1"] = s5
    params_Signal["Midplane"]["s2"] = s6
    params_Signal["Midplane"]["C1"] = C3

    params_Signal["Inplane"]["A1"] = A7
    params_Signal["Inplane"]["A2"] = A8
    params_Signal["Inplane"]["s1"] = s7
    params_Signal["Inplane"]["s2"] = s8
    params_Signal["Inplane"]["C1"] = C4

    
    params_BKG["B"] = B
    params_BKG["v2_t"] = v2_t
    params_BKG["v2_a"] = v2_a
    params_BKG["V3"]   = V3
    params_BKG["v4_t"] = v4_t
    params_BKG["v4_a"] = v4_a

    return Chi2_All(params_Signal, params_BKG) + Chi2_BKG(params_BKG)

def Chi2_All(params_Signal, params_BKG):
    #### FIT TO SIGNAL + BACKGROUND
    total_chi2 = 0
    for key in llaves_plane:
        x = dataALL[key]["x_center"]
        y = dataALL[key]["y"]
        sigma = dataALL[key]["dy"]
        #print x, y
        total_chi2 = total_chi2 + np.sum( np.power( y - Signal( x , params_Signal[key] ) - Background( x , params_BKG, phi_s[key], c[key] ) , 2.0)/ np.power(sigma,2.0) )
    return total_chi2

def Chi2_BKG(params_BKG):
    total_chi2 = 0
    ### BACKGROUND ONLY
    for key in llaves_plane:
        x = dataBKG[key]["x_center"]
        y = dataBKG[key]["y"]
        sigma = dataBKG[key]["dy"]
        length = len(x)
        x = x[:int(length/2)]
        y = y[:int(length/2)]
        sigma = sigma[:int(length/2)]
        #print x, y
        total_chi2 = total_chi2+ np.sum( np.power( y - Background( x , params_BKG, phi_s[key], c[key] )  , 2.0)/ np.power(sigma,2.0))
    return total_chi2

def PerformFitTotal():
    print ' About to start MINUIT'
    sigma_init = 0.07
    limit_sigmaup = 0.5
    m = Minuit(Chi2,
               s3=sigma_init, limit_s3=(0.1,limit_sigmaup), error_s3=0.001,
               s4=sigma_init, limit_s4=(0.1,limit_sigmaup), error_s4=0.001,
               s5=sigma_init, limit_s5=(0.1, limit_sigmaup), error_s5=0.001, 
               s6=sigma_init, limit_s6=(0.1, limit_sigmaup), error_s6=0.001,
               s7=sigma_init, limit_s7=(0.1, limit_sigmaup), error_s7=0.001, 
               s8=sigma_init,  limit_s8=(0.1, limit_sigmaup), error_s8=0.001,
               C2=0.0, fix_C2=True, 
               C3=0.0, fix_C3=True,
               C4=0.0, fix_C4=True,
               #A3 =0.55, limit_A3=(0.1,1.0), error_A3=0.01, 
               #A4=0.05, limit_A4=(0.01, 1.0), error_A4=0.0001,
               #A5 =0.55, limit_A5=(0.1,1.0), error_A5=0.01, 
               #A6=0.05, limit_A6=(0.01, 1.0), error_A6=0.0010,
               #A7 =0.55, limit_A7=(0.1,1.0), error_A7=0.01, 
               #A8=0.05, limit_A8=(0.01, 1.0), error_A8=0.0001,
               v2_t=0.02, limit_v2_t =(0,0.50), error_v2_t=0.001, 
               v2_a=0.02, limit_v2_a =(0,0.50), error_v2_a=0.001,
               v4_t=0.01, limit_v4_t =(0,0.10), error_v4_t=0.001, 
               v4_a=0.01, limit_v4_a =(0,0.10), error_v4_a=0.001,
               V3=0 , limit_V3 = (0, 0.01), error_V3 =0.001)
               #B=2.0 , limit_B = (0.1, 10.0), error_B = 0.01)
    m.migrad()


    fitresult_BKG = {}
    fitresult_BKG["B"]    = m.values["B"]
    fitresult_BKG["v2_t"] = m.values["v2_t"]
    fitresult_BKG["v4_t"] = m.values["v4_t"]
    fitresult_BKG["v2_a"] = m.values["v2_a"]
    fitresult_BKG["v4_a"] = m.values["v4_a"]
    fitresult_BKG["V3"]   = m.values["V3"]
 
    fitresult_Signal = {}

    fitresult_Signal["Outplane"] = {}
    fitresult_Signal["Outplane"]["A1"] = m.values["A3"]
    fitresult_Signal["Outplane"]["A2"] = m.values["A4"]
    fitresult_Signal["Outplane"]["s1"] = m.values["s3"]
    fitresult_Signal["Outplane"]["s2"] = m.values["s4"]
    fitresult_Signal["Outplane"]["C1"] = m.values["C2"]

    fitresult_Signal["Midplane"] = {}
    fitresult_Signal["Midplane"]["A1"] = m.values["A5"]
    fitresult_Signal["Midplane"]["A2"] = m.values["A6"]
    fitresult_Signal["Midplane"]["s1"] = m.values["s5"]
    fitresult_Signal["Midplane"]["s2"] = m.values["s6"]
    fitresult_Signal["Midplane"]["C1"] = m.values["C3"]

    fitresult_Signal["Inplane"] = {}
    fitresult_Signal["Inplane"]["A1"] = m.values["A7"]
    fitresult_Signal["Inplane"]["A2"] = m.values["A8"]
    fitresult_Signal["Inplane"]["s1"] = m.values["s7"]
    fitresult_Signal["Inplane"]["s2"] = m.values["s8"]
    fitresult_Signal["Inplane"]["C1"] = m.values["C4"]


    f, axes = plt.subplots(1,3, sharex=True, sharey=True)


    xfit = np.linspace(-0.5, 0.5, num=150, endpoint=True)
    xfit_all = np.linspace(-0.5, 1.5, num=150, endpoint=True)
    ##FIT RESULT
    for j, key in enumerate(llaves_plane):
        model = Background(xfit, fitresult_BKG, phi_s[key], c[key])
        axes[j].plot(xfit, model,'-', color='dodgerblue') ##plot fit function
        model = TotalPDF(xfit_all, fitresult_Signal[key], fitresult_BKG, phi_s[key], c[key])
        axes[j].plot(xfit_all, model,'-', color='darkorange') ##plot fit function
        x = dataALL[key]["x_center"]
        y = dataALL[key]["y"]
        dy = dataALL[key]["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='o', label='All, \n $|\Delta\eta|<0.5$', color = 'darkorange')
        x = dataBKG[key]["x_center"]#
        y = dataBKG[key]["y"]
        dy = dataBKG[key]["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='o', label='All, \n $1.0<|\Delta\eta|<1.5$', color='dodgerblue')
        axes[j].set_title(key)
        axes[j].set_xlabel(r' $\Delta\phi/\pi$')

        axes[j].locator_params(axis='x', nticks=3)
    axes[0].set_ylabel(r' $1/N_{trig}$ $dN/\Delta\phi$')
    axes[0].legend(loc='best', borderaxespad=0., frameon=False)
    f.subplots_adjust(hspace=0, wspace=0)
    plt.subplots_adjust(hspace=0, wspace=0)
    f.savefig('fitresult_A.png')

    ##SHOW SIGNAL 
    f, axes = plt.subplots(1,3, sharex=True, sharey=True)
    xfit = np.linspace(-0.5, 1.5, num=150, endpoint=True)
    for j, key in enumerate(llaves_plane):
        modelSignal = Signal(xfit, fitresult_Signal[key])
        axes[j].plot(xfit, modelSignal, '-', color='tomato')
        x = dataSIGNAL[key]["x_center"]
        y = dataSIGNAL[key]["y"]
        dy = dataSIGNAL[key]["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='o', label='Signal,\n $|\Delta\eta|<0.5$', color='tomato')
        axes[j].set_title(key)
        axes[j].set_xlabel(r' $\Delta\phi/\pi$')
        s_largeEta = FromTH1toArray(histo["Baseline"]["dphi"]["Signal"]["LargeEta"][key])
        x = s_largeEta["x_center"]
        y = s_largeEta["y"]
        dy = s_largeEta["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='o', label='Signal,\n $1.0<|\Delta\eta|<1.5$', color='lime')
        axes[j].locator_params(axis='x', nticks=4, nbins=4)
        
        
    axes[0].set_ylabel(r' $1/N_{trig}$ $dN/\Delta\phi$')
    axes[0].legend(loc='best', borderaxespad=0., frameon=False)
    plt.show()
    f.subplots_adjust(hspace=0, wspace=0)
    plt.subplots_adjust(hspace=0, wspace=0)
    #plt.tight_layout()
    f.savefig('figresult_B.png')
    ## SHOW BACKGROUND
    f, axes = plt.subplots(1,3, sharex=True, sharey=True)
    xfit = np.linspace(-0.5, 0.5, num=150, endpoint=True)
    for j, key in enumerate(llaves_plane):
        x = dataBKG[key]["x_center"]
        y = dataBKG[key]["y"]
        dy = dataBKG[key]["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='o', label='All \n $1.0<|\Delta\eta|<1.5$', color='crimson')
        modelBKG = Background(xfit, fitresult_BKG, phi_s[key], c[key])
        axes[j].plot(xfit, modelBKG, '-', color='crimson')
        x = dataBKG_True[key]["x_center"]
        y = dataBKG_True[key]["y"]
        dy = dataBKG_True[key]["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='o', label='Bkg \n $1.0<|\Delta\eta|<1.5$', color='lightgrey')
        s_largeEta = FromTH1toArray(histo["Baseline"]["dphi"]["Signal"]["LargeEta"][key])
        x = s_largeEta["x_center"]
        y = s_largeEta["y"]
        dy = s_largeEta["dy"]
        axes[j].errorbar(x,y,yerr=dy, fmt='co', label='Signal, \n $1.0<|\Delta\eta|<1.5$', color='black')
      
        axes[j].set_title(key)
        axes[j].set_xlabel(r' $\Delta\phi/\pi$')
        axes[j].locator_params(axis='x', nticks=3)
        #axes[j].set_ylabel(r' $1/N_{trig}$ $dN/\Delta\phi$')
    axes[0].set_ylabel(r' $1/N_{trig}$ $dN/\Delta\phi$')
    axes[0].legend(loc='best', borderaxespad=0., frameon=False)
    plt.show()
    f.subplots_adjust(hspace=0, wspace=0)
    plt.subplots_adjust(hspace=0, wspace=0)
    f.savefig('figresult_C.png')
   
print 'ola mundo'
#PerformFitToSignal()
PerformFitTotal()
#PerformFitToBKG()
