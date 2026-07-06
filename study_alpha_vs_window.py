#! /usr/bin/python3

import math
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)


plot_folder = '/var/www/html/202602_PDU_cold_test/'

base_folder = '/home/carla/Documents/ADC_data/Legend/20260211_PDU_cold_test_45_7V/run5/'
label = '20260211_45_7V'
ch = 'Ch0'
windows = [500, 1000, 2000, 3000, 4000, 5000, 6000] # in ns
nBins = 500
xMin = 0.
xMax = 2000000.




trees = {}
histos = {}
fits_peak = {}

peak = {}
peak_err = {}
sigma = {}
sigma_err = {}

outfile = ROOT.TFile.Open('%s/%s_alpha_vs_window.root'%(plot_folder,label), 'RECREATE')

for window in windows:
    f = ROOT.TFile.Open(base_folder + 'analyzed_%dns.root'%window)
    trees[window] = f.Get('t1')
    histos[window] = ROOT.TH1F('h_%dns'%window, '', nBins, xMin, xMax)
    trees[window].Draw('%s.charge >> h_%dns'% (ch,window),'','goff')

    spec = ROOT.TSpectrum(2)
    n_pks = spec.Search(histos[window], 10, '',0.01)
    x_pos = spec.GetPositionX()                                                                                                                                                                                                                              
    x_pos = np.array([x_pos[i] for i in range(n_pks)])
    x_pos.sort()
    print(x_pos)
    #fit the last peak with a crystal ball
    if n_pks > 0:
        x_peak = x_pos[-1]
        fit_gaus = ROOT.TF1('fit_gaus_%dns'%window, 'gaus', x_peak-10000, x_peak+10000)
        histos[window].Fit(fit_gaus, 'NSRQ', '', x_peak-100000, x_peak+100000)
        fits_peak[window] = ROOT.TF1('fit_%dns'%window, 'crystalball', 0.75*x_peak, 1.25*x_peak)
        fits_peak[window].SetParameters(fit_gaus.GetParameter(0), fit_gaus.GetParameter(1), fit_gaus.GetParameter(2), 1., 0.5)
        histos[window].Fit(fits_peak[window], 'SRQ+', '', 0.75*x_peak, 1.25*x_peak)
    peak[window] = fits_peak[window].GetParameter(1)
    peak_err[window] = fits_peak[window].GetParError(1)
    sigma[window] = fits_peak[window].GetParameter(2)
    sigma_err[window] = fits_peak[window].GetParError(2)

    c1 = ROOT.TCanvas('c1_%dns'%window, '', 800, 600)
    histos[window].Draw()
    fits_peak[window].Draw('same')
    outfile.cd()
    histos[window].Write()
    c1.SaveAs('%s/%s_%s_alpha_charge_%04dns.png'%(plot_folder,label,ch,window))


c3 = ROOT.TCanvas('c3', '', 800, 600)
hpad = ROOT.gPad.DrawFrame(0., 0., 7000., 1200000.)
hpad.Draw()
hpad.GetXaxis().SetTitle('Integration window [ns]')
hpad.GetXaxis().SetTitleSize(0.05)
hpad.GetXaxis().SetTitleOffset(0.90)
hpad.GetXaxis().SetLabelSize(0.04)
hpad.GetYaxis().SetTitle('#alpha peak [ADC]')
hpad.GetYaxis().SetTitleSize(0.05)
hpad.GetYaxis().SetTitleOffset(1.05)
hpad.GetYaxis().SetLabelSize(0.04)
hpad.GetYaxis().SetTitleColor(ROOT.kRed)
hpad.GetYaxis().SetLabelColor(ROOT.kRed)
ROOT.gPad.SetGridy()
ROOT.gPad.Update()

g1 = ROOT.TGraphErrors(len(windows))
g2 = ROOT.TGraphErrors(len(windows))
for i, window in enumerate(windows):
    g1.SetPoint(i, window, peak[window])
    g1.SetPointError(i, 0, peak_err[window])
    g2.SetPoint(i, window, sigma[window]/peak[window])
    g2.SetPointError(i, 0, sigma_err[window]/peak[window])

g1.SetMarkerStyle(20)
g1.SetMarkerColor(ROOT.kRed)
g1.SetLineColor(ROOT.kRed)
g1.Draw('PL,same')

# asse destro
y1_min = hpad.GetYaxis().GetXmin()
y1_max = hpad.GetYaxis().GetXmax()

y2_min = 0.
y2_max = 0.2

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
axis.SetTitle("peak resolution")
axis.SetTitleColor(ROOT.kBlue)
axis.SetTitleSize(0.05)
axis.SetTitleFont(42)
axis.SetTitleOffset(0.90)
axis.SetLabelSize(0.04)
axis.SetLabelFont(42)
axis.Draw()

c3.SaveAs('%s/%s_%s_alpha_vs_window.png'%(plot_folder,label,ch))

'''
outfile.cd()
g1.Write('gain_vs_window')
g2.Write('SNR_vs_window')
outfile.Close()
'''