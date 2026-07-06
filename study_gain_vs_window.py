#! /usr/bin/python3

import math
import numpy as np
import ROOT
import tdrstyle                                                                                                                                                                        
                                                                                                                                                                                       
#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning






#plot_folder = '/var/www/html/202602_PDU_cold_test/'
#Vbd = 42.09 # in V, extrapoated from LN2 to LAr with 44 mV/K

'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260211_PDU_cold_test_47_2V/run13/'
label = '20260211_47_2V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000, 4000] # in ns
nBins = 1000
xMin = -5000.
xMax = 25000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260211_PDU_cold_test_45_7V/run12/'
label = '20260211_45_7V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260211_PDU_cold_test_45_3V/run12/'
label = '20260211_45_3V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260211_PDU_cold_test_44_7V/run13/'
label = '20260211_44_7V'
ch = 'Ch0'
windows = [500, 750, 1000, 1500, 2000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''
###################################################################################################################
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260213_PDU_cold_test_47_2V/run12/'
label = '20260213_47_2V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000, 4000] # in ns
nBins = 1000
xMin = -5000.
xMax = 25000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260213_PDU_cold_test_45_7V/run12/'
label = '20260213_45_7V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260213_PDU_cold_test_45_3V/run12/'
label = '20260213_45_3V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260213_PDU_cold_test_44_7V/run12/'
label = '20260213_44_7V'
ch = 'Ch0'
windows = [500, 750, 1000, 1500, 2000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''






plot_folder = '/var/www/html/202604_PDU_ModuLYNO_ptp_blu/'
Vbd = 42.09 # in V, extrapoated from LN2 to LAr with 44 mV/K

'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260416_PDU_ModuLYNO_ptp_blu_45_7V/run12/'
label = '20260416_45_7V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260416_PDU_ModuLYNO_ptp_blu_45_83V/run12/'
label = '20260416_45_83V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260416_PDU_ModuLYNO_ptp_blu_47_2V/run12/'
label = '20260416_47_2V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''
'''
base_folder = '/home/carla/Documents/ADC_data/Legend/20260417_PDU_ModuLYNO_ptp_blu_44_7V/run12/'
label = '20260417_44_7V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.
'''

base_folder = '/home/carla/Documents/ADC_data/Legend/20260417_PDU_ModuLYNO_ptp_blu_45_83V/run12/'
label = '20260417_45_83V'
ch = 'Ch0'
windows = [400, 500, 750, 1000, 1500, 2000, 3000] # in ns
nBins = 1000
xMin = -5000.
xMax = 20000.


trees = {}
histos = {}
fits_sphe = {}

gain = {}
gain_err = {}
SNR = {}

outfile = ROOT.TFile.Open('%s/%s_gain_vs_window.root'%(plot_folder,label), 'RECREATE')

Vbias = float(label.split('_')[1]+'.'+label.split('_')[2].replace('V',''))
Vov = Vbias - Vbd - 1.

for window in windows:
    f = ROOT.TFile.Open(base_folder + 'analyzed_%dns.root'%window)
    trees[window] = f.Get('t1')
    histos[window] = ROOT.TH1F('h_%dns'%window, '', nBins, xMin, xMax)
    trees[window].Draw('%s.charge >> h_%dns'% (ch,window),'','goff')

    spec = ROOT.TSpectrum(5)
    n_pks = spec.Search(histos[window], 10, '',0.01)
    x_pos = spec.GetPositionX()                                                                                                                                                                                                                              
    x_pos = np.array([x_pos[i] for i in range(n_pks)])
    x_pos.sort()
    print(x_pos)
    for i in range(n_pks):
        delta = x_pos[1] - x_pos[0]
        fits_sphe[(window,i)] = ROOT.TF1('fit_%dns_%d'%(window,i), 'gaus', x_pos[i]-0.25*delta, x_pos[i]+0.25*delta)
        histos[window].Fit(fits_sphe[(window,i)], 'SRQ+', '', x_pos[i]-0.25*delta, x_pos[i]+0.25*delta)
    gain[window] = fits_sphe[(window,2)].GetParameter(1) - fits_sphe[(window,1)].GetParameter(1)
    gain_err[window] = math.sqrt(fits_sphe[(window,2)].GetParError(1)**2 + fits_sphe[(window,1)].GetParError(1)**2)
    SNR[window] = ( fits_sphe[(window,1)].GetParameter(1) ) / ( fits_sphe[(window,1)].GetParameter(2) )
    print('Window: %d ns, Gain: %.2f +/- %.2f, SNR: %.2f'%(window, gain[window], gain_err[window], SNR[window]))

    c1 = ROOT.TCanvas('c1_%dns'%window, '', 800, 600)
    histos[window].SetTitle(';Charge [ADC];Counts')
    histos[window].Draw()
    outfile.cd()
    histos[window].Write()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.05)
    latex.SetTextFont(42)
    latex.DrawLatex(0.65, 0.84, 'V_{OV} = %.2f V'%(Vov))
    c1.SaveAs('%s/%s_%s_charge_%04dns.png'%(plot_folder,label,ch,window))

#plot gain and SNR vs window on the same graph (two different y axes, one for gain and one for SNR)
c3 = ROOT.TCanvas('c3', '', 800, 600)
hpad = ROOT.gPad.DrawFrame(0., 0., 5000., 6000.)
hpad.Draw()
hpad.GetXaxis().SetTitle('Integration window [ns]')
hpad.GetXaxis().SetTitleSize(0.05)
hpad.GetXaxis().SetTitleOffset(0.90)
hpad.GetXaxis().SetLabelSize(0.04)
hpad.GetYaxis().SetTitle('Gain [a.u.]')
hpad.GetYaxis().SetTitleSize(0.05)
hpad.GetYaxis().SetTitleOffset(1.05)
hpad.GetYaxis().SetLabelSize(0.04)
hpad.GetYaxis().SetTitleColor(ROOT.kRed)
hpad.GetYaxis().SetLabelColor(ROOT.kRed)
ROOT.gPad.SetGridy()
ROOT.gPad.Update()

g1 = ROOT.TGraphErrors(len(windows))
g2 = ROOT.TGraph(len(windows))
for i, window in enumerate(windows):
    g1.SetPoint(i, window, gain[window])
    g1.SetPointError(i, 0, gain_err[window])
    g2.SetPoint(i, window, SNR[window])

g1.SetMarkerStyle(20)
g1.SetMarkerColor(ROOT.kRed)
g1.SetLineColor(ROOT.kRed)
g1.Draw('PL,same')

# asse destro
y1_min = hpad.GetYaxis().GetXmin()
y1_max = hpad.GetYaxis().GetXmax()

y2_min = 0.  #yvals.min() - 0.1 * (yvals.max() - yvals.min())
y2_max = 15. #yvals.max() + 0.1 * (yvals.max() - yvals.min())

scale = (y1_max - y1_min) / ( y2_max - y2_min )

n = g2.GetN()
g2s = ROOT.TGraph(n)
yvals = np.array([g2.GetPointY(i) for i in range(n)])
for i in range(n):
    x = g2.GetPointX(i)
    y = g2.GetPointY(i)
    yy = y1_min + (y - y2_min) * scale
    g2s.SetPoint(i, x, yy)

g2s.SetMarkerStyle(20)
g2s.SetMarkerColor(ROOT.kBlue)
g2s.SetLineColor(ROOT.kBlue)
g2s.Draw("PL SAME")
ROOT.gPad.Update()

axis = ROOT.TGaxis(
    ROOT.gPad.GetUxmax(), y1_min,
    ROOT.gPad.GetUxmax(), y1_max,
    y2_min, y2_max,
    510, "+L"
)
axis.SetLineColor(ROOT.kBlue)
axis.SetLabelColor(ROOT.kBlue)
axis.SetTitle("S / N")
axis.SetTitleColor(ROOT.kBlue)
axis.SetTitleSize(0.05)
axis.SetTitleFont(42)
axis.SetTitleOffset(0.80)
axis.SetLabelSize(0.04)
axis.SetLabelFont(42)
axis.Draw()

c3.SaveAs('%s/%s_%s_gain_and_SNR_vs_window.png'%(plot_folder,label,ch))

outfile.cd()
g1.Write('gain_vs_window')
g2.Write('SNR_vs_window')
outfile.Close()
