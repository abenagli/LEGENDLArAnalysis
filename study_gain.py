#! /usr/bin/python3

import os
import math
import numpy as np
import ROOT
import tdrstyle                                                                                                                                                                        
                                                                                                                                                                                       
#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

#Vbd = 42.09 # in V, extrapoated from LN2 to LAr with 44 mV/K
#plot_folder = '/var/www/html/202602_PDU_cold_test/'
#labels = ['20260211', '20260213']
#Vbias = ['47_2V', '45_7V', '45_3V', '44_7V']

Vbd = 42.09 # in V, extrapoated from LN2 to LAr with 44 mV/K
plot_folder = '/var/www/html/202604_PDU_ModuLYNO_ptp_blu/'
labels = ['20260416', '20260417']
Vbias = ['44_7V', '45_7V', '45_83V', '47_2V']

infiles = {}
graphs_gain = {}
graphs_SNR = {}

g_gain_vs_Vov = {}
g_SNR_vs_Vov = {}

for label in labels:
    outfile = ROOT.TFile.Open('%s/%s_gain.root'%(plot_folder,label), 'RECREATE')

    g_gain_vs_Vov[label] = ROOT.TGraphErrors()
    g_SNR_vs_Vov[label] = ROOT.TGraphErrors()

    for V in Vbias:
        if not os.path.exists('%s/%s_%s_gain_vs_window.root'%(plot_folder,label,V)):
            print("File %s/%s_%s_gain_vs_window.root does not exist. Skipping."%(plot_folder,label,V))
            continue
        infiles[label+'_'+V] = ROOT.TFile.Open('%s/%s_%s_gain_vs_window.root'%(plot_folder,label,V))
        graphs_gain[label+'_'+V] = infiles[label+'_'+V].Get('gain_vs_window')
        graphs_SNR[label+'_'+V] = infiles[label+'_'+V].Get('SNR_vs_window')

    for V in Vbias:
        if label+'_'+V not in graphs_gain:
            print("No graph for %s at %s. Skipping."%(label, V))
            continue
        f_corr = graphs_gain[label+'_'+V].Eval(3000.) / graphs_gain[label+'_'+V].Eval(1000.)
        print(f_corr)

        Vov = float((V[:-1]).replace('_', '.')) - Vbd - 1.
        gain = graphs_gain[label+'_'+V].Eval(1000.)
        gain_err = graphs_gain[label+'_'+V].GetErrorY(3)
        SNR = graphs_SNR[label+'_'+V].Eval(1000.)
        g_gain_vs_Vov[label].SetPoint(g_gain_vs_Vov[label].GetN(), Vov, gain*f_corr)
        g_gain_vs_Vov[label].SetPointError(g_gain_vs_Vov[label].GetN()-1, 0, gain_err*f_corr)
        g_SNR_vs_Vov[label].SetPoint(g_SNR_vs_Vov[label].GetN(), Vov, SNR)
        g_SNR_vs_Vov[label].SetPointError(g_SNR_vs_Vov[label].GetN()-1, 0, 0)
    
    outfile.cd()
    g_gain_vs_Vov[label].Write('gain_vs_Vov_%s'%label)
    g_SNR_vs_Vov[label].Write('SNR_vs_Vov_%s'%label)

outfile.Close()


c1 = ROOT.TCanvas('c_gain_SNR_vs_Vov', 'c_gain_SNR_vs_Vov', 800, 700)
hpad = ROOT.gPad.DrawFrame(0., 0., 5., 6000.)
hpad.Draw()
hpad.SetTitle(';V_{ov} [V];Gain [a.u.]')
hpad.GetYaxis().SetTitleColor(ROOT.kRed)
hpad.GetYaxis().SetLabelColor(ROOT.kRed)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetTicky(0)
ROOT.gPad.Update()

it = 0
for label in labels:
    if it == 0: g_gain_vs_Vov[label].SetMarkerStyle(20)
    else: g_gain_vs_Vov[label].SetMarkerStyle(24)
    g_gain_vs_Vov[label].SetMarkerSize(1.3)
    g_gain_vs_Vov[label].SetMarkerColor(ROOT.kRed)
    g_gain_vs_Vov[label].SetLineColor(ROOT.kRed)
    if it == 1: g_gain_vs_Vov[label].SetLineStyle(7)
    g_gain_vs_Vov[label].Draw('P,same')
    
    fit_gain = ROOT.TF1('fit_gain_%s'%(label), '[0]*(x+[1])', 0., 5.)
    fit_gain.SetLineColor(ROOT.kRed)
    fit_gain.SetLineWidth(2)
    if it == 0: fit_gain.SetLineStyle(1)
    else: fit_gain.SetLineStyle(7)
    g_gain_vs_Vov[label].Fit(fit_gain, 'QNRS+')
    fit_gain.Draw('same')
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextColor(ROOT.kRed)
    latex.DrawLatex(0.17, 0.80+0.05*it, 'Gain = %.1f #times ( V_{ov} + %.2f V)'%(fit_gain.GetParameter(0), fit_gain.GetParameter(1)))
    it += 1

# asse destro
y1_min = hpad.GetYaxis().GetXmin()
y1_max = hpad.GetYaxis().GetXmax()

y2_min = 0.
y2_max = 20.

scale = (y1_max - y1_min) / ( y2_max - y2_min )

g2s = {}
it = 0
for label in labels:
    n = g_SNR_vs_Vov[label].GetN()
    g2s[label] = ROOT.TGraph(n)
    for i in range(n):
        x = g_SNR_vs_Vov[label].GetPointX(i)
        y = g_SNR_vs_Vov[label].GetPointY(i)
        yy = y1_min + (y - y2_min) * scale
        g2s[label].SetPoint(i, x, yy)

    if it == 0: g2s[label].SetMarkerStyle(20)
    else: g2s[label].SetMarkerStyle(24)
    g2s[label].SetMarkerSize(1.3)
    g2s[label].SetMarkerColor(ROOT.kBlue)
    g2s[label].SetLineColor(ROOT.kBlue)
    g2s[label].SetLineWidth(2)
    if it == 0: g2s[label].SetLineStyle(1)
    else: g2s[label].SetLineStyle(7)
    g2s[label].Draw("PL SAME")
    ROOT.gPad.Update()
    it += 1

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
axis.SetLabelSize(0.05)
axis.SetLabelFont(42)
axis.Draw()

legend = ROOT.TLegend(0.35, 0.18, 0.85, 0.30)
legend.SetBorderSize(0)
legend.SetFillStyle(1001)
legend.SetTextSize(0.05)
legend.SetTextFont(42)
#legend.AddEntry(g_gain_vs_Vov[labels[0]], '10 SiPMs ganged', 'PL')
#legend.AddEntry(g_gain_vs_Vov[labels[1]], '5 SiPMs ganged.', 'PL')
legend.AddEntry(g_gain_vs_Vov[labels[0]], '2026-04-16', 'PL')
legend.AddEntry(g_gain_vs_Vov[labels[1]], '2026-04-17', 'PL')
legend.Draw()

c1.SaveAs('%s/gain_vs_Vov.png'%(plot_folder))
