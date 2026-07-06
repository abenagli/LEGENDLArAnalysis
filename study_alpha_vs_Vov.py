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


Vbd = 42.09 # in V, extrapoated from LN2 to LAr with 44 mV/K

plot_folder = '/var/www/html/202602_PDU_cold_test/'

labels = ['20260211_45_7V', '20260211_45_3V', '20260211_44_7V']
labels2 = ['20260213_45_7V', '20260213_45_3V']
base_folder = '/home/carla/Documents/ADC_data/Legend/20260213_PDU_cold_test_45_7V/'

g_phe_vs_Vov = ROOT.TGraphErrors()
g_phe2_vs_Vov = ROOT.TGraphErrors()

for label in labels:
    infile = ROOT.TFile.Open('/var/www/html/202602_PDU_cold_test/%s_alpha.root'%label, 'READ')
    g_sum = infile.Get('g_sum')
    fit_pol0 = ROOT.TF1('fit_pol0', 'pol0', 7., 41.)
    g_sum.Fit(fit_pol0, 'RQ')

    Vbias = float(label.split('_')[1]+'.'+label.split('_')[2].replace('V',''))
    Vov = Vbias - Vbd - 1.
    g_phe_vs_Vov.SetPoint(g_phe_vs_Vov.GetN(), Vov, fit_pol0.GetParameter(0))
    g_phe_vs_Vov.SetPointError(g_phe_vs_Vov.GetN()-1, 0., fit_pol0.GetParError(0))

for label in labels2:
    infile = ROOT.TFile.Open('/var/www/html/202602_PDU_cold_test/%s_alpha.root'%label, 'READ')
    g_sum = infile.Get('g_sum')
    fit_pol0 = ROOT.TF1('fit_pol0', 'pol0', 7., 41.)
    g_sum.Fit(fit_pol0, 'RQ')

    Vbias = float(label.split('_')[1]+'.'+label.split('_')[2].replace('V',''))
    Vov = Vbias - Vbd - 1.
    g_phe2_vs_Vov.SetPoint(g_phe2_vs_Vov.GetN(), Vov, fit_pol0.GetParameter(0))
    g_phe2_vs_Vov.SetPointError(g_phe2_vs_Vov.GetN()-1, 0., fit_pol0.GetParError(0))

c1 = ROOT.TCanvas('c1', '', 800, 600)
hpad = ROOT.gPad.DrawFrame(0., 0., 5., 800.)
hpad.GetXaxis().SetTitle('V_{ov} [V]')
hpad.GetYaxis().SetTitle('#alpha peak [N_{ph.e.}]')
g_phe_vs_Vov.SetMarkerStyle(20)
g_phe_vs_Vov.SetMarkerColor(ROOT.kBlack)
g_phe_vs_Vov.SetLineColor(ROOT.kBlack)
g_phe2_vs_Vov.SetMarkerStyle(24)
g_phe2_vs_Vov.SetMarkerColor(ROOT.kBlack)
g_phe2_vs_Vov.SetLineColor(ROOT.kBlack)
g_phe2_vs_Vov.SetLineStyle(7)
g_phe_vs_Vov.Draw('PL,same')
g_phe2_vs_Vov.Draw('PL,same')

legend = ROOT.TLegend(0.38, 0.18, 0.85, 0.30)
legend.SetBorderSize(0)
legend.SetFillStyle(1001)
legend.SetTextSize(0.05)
legend.SetTextFont(42)
legend.AddEntry(g_phe_vs_Vov, '10 SiPMs connected', 'PL')
legend.AddEntry(g_phe2_vs_Vov, '5 SiPMs connected', 'PL')
legend.Draw()

c1.SaveAs('%s/phe_vs_Vov.png'%(plot_folder))
