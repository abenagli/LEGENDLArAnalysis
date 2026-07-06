#! /usr/bin/python3

import math
import numpy as np
import ROOT
import tdrstyle
import os

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


Vbd = 42.09 # in V, extrapoated from LN2 to LAr with 44 mV/K

# configurable fprompt cut
fprompt_cut = 0.48

# default charge_1pe if no gain file is provided
charge_1pe_default = 1.0

positions = {}
positions['run1'] = 4.
positions['run2'] = 8.
positions['run3'] = 12.
positions['run41'] = 20.
positions['run51'] = 24.
positions['run61'] = 28.
positions['run7'] = 36.
positions['run8'] = 40.
positions['run9'] = 44.



'''
plot_folder = '/var/www/html/202602_PDU_cold_test/'

base_folder = '/home/carla/Documents/ADC_data/Legend/20260211_PDU_cold_test_45_7V/'
label = '20260211_45_7V'
file_gain_name = '/var/www/html/202602_PDU_cold_test/20260211_gain.root'
#base_folder = '/home/carla/Documents/ADC_data/Legend/20260211_PDU_cold_test_45_3V/'
#label = '20260211_45_3V'
#file_gain_name = '/var/www/html/202602_PDU_cold_test/20260211_gain.root'
#base_folder = '/home/carla/Documents/ADC_data/Legend/20260211_PDU_cold_test_44_7V/'
#label = '20260211_44_7V'
#file_gain_name = '/var/www/html/202602_PDU_cold_test/20260211_gain.root'

#base_folder = '/home/carla/Documents/ADC_data/Legend/20260213_PDU_cold_test_45_7V/'
#label = '20260213_45_7V'
#file_gain_name = '/var/www/html/202602_PDU_cold_test/20260213_gain.root'
#base_folder = '/home/carla/Documents/ADC_data/Legend/20260213_PDU_cold_test_45_3V/'
#label = '20260213_45_3V'
#file_gain_name = '/var/www/html/202602_PDU_cold_test/20260213_gain.root'
'''



#plot_folder = '/var/www/html/202604_PDU_ModuLYNO_ptp_blu/'

#base_folder = '/home/carla/Documents/ADC_data/Legend/20260416_PDU_ModuLYNO_ptp_blu_45_7V/'
#label = '20260416_45_7V'
#file_gain_name = '/var/www/html/202604_PDU_ModuLYNO_ptp_blu/20260416_gain.root'
#base_folder = '/home/carla/Documents/ADC_data/Legend/20260416_PDU_ModuLYNO_ptp_blu_45_83V/'
#label = '20260416_45_83V'
#file_gain_name = '/var/www/html/202604_PDU_ModuLYNO_ptp_blu/20260416_gain.root'
#base_folder = '/home/carla/Documents/ADC_data/Legend/20260416_PDU_ModuLYNO_ptp_blu_47_2V/'
#label = '20260416_47_2V'
#file_gain_name = '/var/www/html/202604_PDU_ModuLYNO_ptp_blu/20260416_gain.root'



plot_folder = '/var/www/html/20260605_PDU_ModulYNO_blu-verde_noWindow_laser'

base_folder = '/home/carla/Documents/ADC_data/Legend/20260605_PDU_ModulYNO_blu-verde_noWindow_laser_55_0V/'
label = '20260605_55_0V'
file_gain_name = None



runs = ['run1', 'run2', 'run3', 'run41', 'run51', 'run61', 'run7', 'run8', 'run9']
nBins = 500
xMin = 0.
xMax = 3000000.




trees = {}
histos = {}
fits = {}

peaks = {}
peaks_err = {}
sigmas = {}
sigmas_err = {}

outfile = ROOT.TFile.Open('%s/%s_alpha.root'%(plot_folder,label), 'RECREATE')

Vbias = float(label.split('_')[1]+'.'+label.split('_')[2].replace('V',''))
Vov = Vbias - Vbd - 1.

use_gain_file = False
charge_1pe = charge_1pe_default

if file_gain_name is not None and os.path.exists(file_gain_name):
    file_gain = ROOT.TFile.Open(file_gain_name, 'READ')
    if file_gain and not file_gain.IsZombie():
        gain_name = 'gain_vs_Vov_%s'%(label.split('_')[0])
        gain_obj = file_gain.Get(gain_name)
        if gain_obj:
            charge_1pe = gain_obj.Eval(Vov)
            use_gain_file = True
        else:
            print("Gain object %s not found. Using default charge_1pe = %.3f."%(gain_name, charge_1pe_default))
    else:
        print("Could not open gain file %s. Using default charge_1pe = %.3f."%(file_gain_name, charge_1pe_default))
else:
    print("No gain file provided. Using default charge_1pe = %.3f."%charge_1pe_default)

y_axis_label = '#alpha peak [N_{ph.e.}]' if use_gain_file else '#alpha peak [ADC]'

for run in runs:
    #check if file exists
    if not os.path.exists(base_folder + run + '/analyzed.root'):
        print("File %s does not exist. Skipping."%(base_folder + run + '/analyzed.root'))
        continue
    f = ROOT.TFile.Open(base_folder + run + '/analyzed.root')
    trees[run] = f.Get('t1')
    histos[(run,'ch0')] = ROOT.TH1F('h_ch0_%s'%run, '', nBins, xMin, xMax)
    histos[(run,'ch1')] = ROOT.TH1F('h_ch1_%s'%run, '', nBins, xMin, xMax)
    histos[(run,'sum')] = ROOT.TH1F('h_sum_%s'%run, '', nBins, xMin, xMax)

    cut_ch0 = 'Ch0.fprompt>%.6g'%fprompt_cut
    cut_ch1 = 'Ch1.fprompt>%.6g'%fprompt_cut
    cut_sum = 'Ch0.fprompt>%.6g && Ch1.fprompt>%.6g'%(fprompt_cut, fprompt_cut)

    trees[run].Draw('Ch0.charge >> h_ch0_%s'% run, cut_ch0, 'goff')
    trees[run].Draw('Ch1.charge >> h_ch1_%s'% run, cut_ch1, 'goff')
    trees[run].Draw('(Ch0.charge + Ch1.charge) >> h_sum_%s'% run, cut_sum, 'goff')

    for ch in ['ch0', 'ch1', 'sum']:
        spec = ROOT.TSpectrum(2)
        n_pks = spec.Search(histos[(run,ch)], 10, '',0.01)
        x_pos = spec.GetPositionX()
        x_pos = np.array([x_pos[i] for i in range(n_pks)])
        x_pos.sort()
        print(x_pos)
        if n_pks > 0:
            x_peak = x_pos[-1]
            fit_gaus = ROOT.TF1('fit_gaus_%s'%run, 'gaus', x_peak-10000, x_peak+10000)
            histos[(run,ch)].Fit(fit_gaus, 'NSRQ', '', x_peak-100000, x_peak+100000)
            fits[(run,ch)] = ROOT.TF1('fit_%s_%s'%(ch,run), 'crystalball', 0.75*x_peak, 1.25*x_peak)
            fits[(run,ch)].SetParameters(fit_gaus.GetParameter(0), fit_gaus.GetParameter(1), fit_gaus.GetParameter(2), 1., 0.5)
            histos[(run,ch)].Fit(fits[(run,ch)], 'SRQ', '', 0.75*x_peak, 1.25*x_peak)
        peaks[(run,ch)] = fits[(run,ch)].GetParameter(1)
        peaks_err[(run,ch)] = fits[(run,ch)].GetParError(1)
        sigmas[(run,ch)] = fits[(run,ch)].GetParameter(2)
        sigmas_err[(run,ch)] = fits[(run,ch)].GetParError(2)


    c1 = ROOT.TCanvas('c1_%s'%run, '', 800, 600)
    histos[(run,'sum')].SetMaximum(1.2*max(histos[(run,'sum')].GetMaximum(), histos[(run,'ch0')].GetMaximum(), histos[(run,'ch1')].GetMaximum()))
    histos[(run,'sum')].SetLineColor(ROOT.kBlack)
    histos[(run,'sum')].SetLineWidth(2)
    histos[(run,'ch0')].SetLineColor(ROOT.kRed)
    histos[(run,'ch0')].SetLineWidth(2)
    histos[(run,'ch1')].SetLineColor(ROOT.kBlue)
    histos[(run,'ch1')].SetLineWidth(2)
    histos[(run,'sum')].SetTitle(';Charge [ADC];Counts')
    histos[(run,'sum')].Draw()
    histos[(run,'ch0')].Draw('same')
    histos[(run,'ch1')].Draw('same')

    fits[(run,'sum')].SetLineColor(ROOT.kBlack)
    fits[(run,'sum')].SetLineWidth(2)
    fits[(run,'ch0')].SetLineColor(ROOT.kRed)
    fits[(run,'ch0')].SetLineWidth(2)
    fits[(run,'ch1')].SetLineColor(ROOT.kBlue)
    fits[(run,'ch1')].SetLineWidth(2)
    fits[(run,'ch0')].Draw('same')
    fits[(run,'ch1')].Draw('same')
    fits[(run,'sum')].Draw('same')

    leg = ROOT.TLegend(0.60, 0.62, 0.88, 0.80)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.AddEntry(histos[(run,'ch0')], 'top', 'l')
    leg.AddEntry(histos[(run,'ch1')], 'bottom', 'l')
    leg.AddEntry(histos[(run,'sum')], 'sum', 'l')
    leg.Draw('same')

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.05)
    latex.SetTextFont(42)
    latex.DrawLatex(0.65, 0.84, 'V_{OV} = %.2f V'%(Vov))
    c1.SaveAs('%s/%s_alpha_charge_%s.png'%(plot_folder,label,run))

    outfile.cd()
    histos[(run,'ch0')].Write()
    histos[(run,'ch1')].Write()
    histos[(run,'sum')].Write()
    fits[(run,'ch0')].Write()
    fits[(run,'ch1')].Write()
    fits[(run,'sum')].Write()


# draw peak vs position
c2 = ROOT.TCanvas('c2', '', 800, 600)
hpad = ROOT.gPad.DrawFrame(0., 0., 50., 1500000.)
hpad.Draw()
ROOT.gPad.SetGridy()
hpad.GetXaxis().SetTitle('Position [mm]')
hpad.GetYaxis().SetTitle(y_axis_label)

g_ch0 = ROOT.TGraphErrors()
g_ch1 = ROOT.TGraphErrors()
g_sum = ROOT.TGraphErrors()
for i, run in enumerate(runs):
    if (run,'ch0') not in peaks or (run,'ch1') not in peaks or (run,'sum') not in peaks:
        print("No peak found for %s. Skipping point."%run)
        continue
    print(g_ch0.GetN(), positions[run], peaks[(run,'ch0')]/charge_1pe, peaks_err[(run,'ch0')]/charge_1pe)
    g_ch0.SetPoint(g_ch0.GetN(), positions[run], peaks[(run,'ch0')]/charge_1pe)
    g_ch0.SetPointError(g_ch0.GetN()-1, 0, peaks_err[(run,'ch0')]/charge_1pe)
    g_ch1.SetPoint(g_ch1.GetN(), positions[run], peaks[(run,'ch1')]/charge_1pe)
    g_ch1.SetPointError(g_ch1.GetN()-1, 0, peaks_err[(run,'ch1')]/charge_1pe)
    g_sum.SetPoint(g_sum.GetN(), positions[run], peaks[(run,'sum')]/charge_1pe)
    g_sum.SetPointError(g_sum.GetN()-1, 0, peaks_err[(run,'sum')]/charge_1pe)

g_sum.SetMarkerStyle(20)
g_sum.SetMarkerColor(ROOT.kBlack)
g_sum.SetLineColor(ROOT.kBlack)
g_ch0.SetMarkerStyle(20)
g_ch0.SetMarkerColor(ROOT.kRed)
g_ch0.SetLineColor(ROOT.kRed)
g_ch1.SetMarkerStyle(20)
g_ch1.SetMarkerColor(ROOT.kBlue)
g_ch1.SetLineColor(ROOT.kBlue)

g_ch0.Draw('PL,same')
g_ch1.Draw('PL,same')
g_sum.Draw('PL,same')

leg2 = ROOT.TLegend(0.60, 0.62, 0.88, 0.80)
leg2.SetBorderSize(0)
leg2.SetFillStyle(0)
leg2.SetTextFont(42)
leg2.SetTextSize(0.04)
leg2.AddEntry(g_ch0, 'top', 'pl')
leg2.AddEntry(g_ch1, 'bottom', 'pl')
leg2.AddEntry(g_sum, 'sum', 'pl')
leg2.Draw('same')

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.05)
latex.SetTextFont(42)
latex.DrawLatex(0.65, 0.84, 'V_{OV} = %.2f V'%(Vov))
c2.SaveAs('%s/%s_alpha_peak_vs_position.png'%(plot_folder,label))

outfile.cd()
g_ch0.Write('g_ch0')
g_ch1.Write('g_ch1')
g_sum.Write('g_sum')
outfile.Close()
