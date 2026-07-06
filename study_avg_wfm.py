#!/usr/bin/env python3
import os
import numpy as np
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)


'''
label = '20251117_LYno_pTPonPMMA'
run = 'run11_45V00_600ADC_Ch0'
plot_folder = '/var/www/html/202602_PDU_cold_test/'
charge_min = 0.
charge_max = 7e6
fprompt_min = 0.
fprompt_max = 0.4
'''



label = '20260211_PDU_cold_test_44_7V'
run = 'run5'
plot_folder = '/var/www/html/202602_PDU_cold_test/'
charge_min = 0.
charge_max = 7e6
fprompt_min = 0.4
fprompt_max = 1.


'''
label = '20260211_PDU_cold_test_44_7V'
run = 'run11'
plot_folder = '/var/www/html/202602_PDU_cold_test/'
charge_min = 0.
charge_max = 7e6
fprompt_min = 0.
fprompt_max = 0.4
'''

'''
label = '20260211_PDU_cold_test_44_7V'
run = 'run13'
plot_folder = '/var/www/html/202602_PDU_cold_test/'
charge_min = 6200.
charge_max = 7200
fprompt_min = 0.
fprompt_max = 1.
'''

bin_width = 4. # in ns
tree_name = "t1"
w_branch = "Ch0.wvf"
q_branch = "Ch0.charge"
fp_branch = "Ch0.fprompt"

base_folder = '/home/carla/Documents/ADC_data/Legend/'

npts = 5000
NBASE = 950
PULSE_SIGN = +1


def t50_crossing(w: np.ndarray, nbase: int = 100, pulse_sign: int = +1):
    """
    Ritorna t50 (float, sample frazionario) usando crossing a metà ampiezza sul fronte di salita.
    Se non trova, ritorna None.
    """
    if w.size == 0:
        return None

    base = float(np.mean(w[:nbase])) if w.size >= nbase else float(np.mean(w))

    if pulse_sign > 0:
        # pulse positiva: usa massimo
        peak_idx = int(np.nanargmax(w))
        peak = float(w[peak_idx])
        amp = peak - base
        if not np.isfinite(amp) or amp <= 0:
            return None
        half = base + 0.5 * amp

        # cerca crossing da sotto a sopra prima del picco
        # (fronte di salita)
        start = 0
        stop = peak_idx if peak_idx > 1 else w.size - 1
        for i in range(start, stop):
            y0, y1 = float(w[i]), float(w[i+1])
            if y0 < half <= y1:
                # interpolazione lineare
                denom = (y1 - y0)
                if denom == 0:
                    return float(i)
                frac = (half - y0) / denom
                return float(i) + float(frac)
        return None
    else:
        # pulse negativa: usa minimo
        peak_idx = int(np.nanargmin(w))
        peak = float(w[peak_idx])
        amp = base - peak
        if not np.isfinite(amp) or amp <= 0:
            return None
        half = base - 0.5 * amp

        # crossing da sopra a sotto prima del minimo
        start = 0
        stop = peak_idx if peak_idx > 1 else w.size - 1
        for i in range(start, stop):
            y0, y1 = float(w[i]), float(w[i+1])
            if y0 > half >= y1:
                denom = (y1 - y0)
                if denom == 0:
                    return float(i)
                frac = (half - y0) / denom
                return float(i) + float(frac)
        return None


def shift_waveform(w: np.ndarray, shift: float):
    """
    Applica shift orizzontale frazionario (in sample).
    shift > 0: sposta a destra (ritarda).
    """
    x = np.arange(w.size, dtype=np.float64)
    xq = x - shift
    # left/right -> NaN così non inquina la media
    return np.interp(xq, x, w, left=np.nan, right=np.nan)




# ---- I/O ROOT
infile = ROOT.TFile.Open(os.path.join(base_folder, label, run, 'analyzed.root'))
if not infile or infile.IsZombie():
    raise RuntimeError("File non apribile")

t = infile.Get(tree_name)
if not t:
    raise RuntimeError(f"Tree {tree_name} non trovato")

t.SetBranchStatus("*", 0)
t.SetBranchStatus(q_branch, 1)
t.SetBranchStatus(w_branch, 1)
t.SetBranchStatus(fp_branch, 1)

q_buf = array('d', [0.0])
fp_buf = array('d', [0.0])
w_buf = np.zeros(npts, dtype=np.float64)

t.SetBranchAddress(q_branch, q_buf)
t.SetBranchAddress(fp_branch, fp_buf)
t.SetBranchAddress(w_branch, w_buf)

n_tot = t.GetEntries()

h_base = ROOT.TH1F('h_base', '', 10000, -100., 100.)
h_charge_mid = ROOT.TH1F('h_charge_mid', '', 1000, -10000., 30000.)
h_charge_late = ROOT.TH1F('h_charge_late', '', 1000, -10000., 30000.)

# ---- PASS 1: calcola t50 per tutti gli eventi selezionati
t50_list = []

for i in range(n_tot):
    t.GetEntry(i)

    base = float(np.mean(w_buf[:NBASE])) if w_buf.size >= NBASE else float(np.mean(w_buf))
    h_base.Fill(base)
    if abs(base) > 0.1:
            continue

    charge_mid = float(np.sum(w_buf[1500:4150]))
    h_charge_mid.Fill(charge_mid)
    #if charge_mid > 4000.:
    #    continue
    
    charge_late = float(np.sum(w_buf[4150:4999]))
    h_charge_late.Fill(charge_late)
    #if charge_late > 1000:
    #    continue
    
    q = float(q_buf[0])
    fp = float(fp_buf[0])

    if not np.isfinite(q) or q == 0.0:
        continue
    if q < charge_min or q > charge_max:
        continue
    if fp < fprompt_min or fp > fprompt_max:
        continue

    # t50 calcolato sulla waveform "raw" (normalizzare non cambia il tempo)
    tt = t50_crossing(w_buf, nbase=NBASE, pulse_sign=PULSE_SIGN)
    if tt is None or not np.isfinite(tt):
        continue
    t50_list.append(tt)

if len(t50_list) == 0:
    raise RuntimeError("Nessun t50 valido trovato (controlla segno pulse / baseline / tagli).")

t_ref = float(np.median(np.asarray(t50_list, dtype=np.float64)))
print(f"t_ref (median t50) = {t_ref:.3f} samples  (N={len(t50_list)})")


# ---- PASS 2: accumula media di (wvf/q) shiftata
sum_bins = np.zeros(npts, dtype=np.float64)
counts   = np.zeros(npts, dtype=np.int64)
n_used = 0

for i in range(n_tot):
    t.GetEntry(i)

    base = float(np.mean(w_buf[:NBASE])) if w_buf.size >= NBASE else float(np.mean(w_buf))
    h_base.Fill(base)
    if abs(base) > 0.1:
            continue
    
    charge_mid = float(np.sum(w_buf[1500:4150]))
    h_charge_mid.Fill(charge_mid)
    #if charge_mid > 4000.:
    #    continue
        
    charge_late = float(np.sum(w_buf[4150:4999]))
    h_charge_late.Fill(charge_late)
    #if charge_late > 20000:
    #    continue

    q = float(q_buf[0])
    fp = float(fp_buf[0])

    if not np.isfinite(q) or q == 0.0:
        continue
    if q < charge_min or q > charge_max:
        continue
    if fp < fprompt_min or fp > fprompt_max:
        continue

    tt = t50_crossing(w_buf, nbase=NBASE, pulse_sign=PULSE_SIGN)
    if tt is None or not np.isfinite(tt):
        continue

    shift = t_ref - tt

    w_norm = w_buf / q
    w_shift = shift_waveform(w_norm, shift)

    valid = np.isfinite(w_shift)
    sum_bins[valid] += w_shift[valid]
    counts[valid]   += 1
    n_used += 1


print(f"Eventi totali: {n_tot}, usati: {n_used}")
if n_used == 0:
    raise RuntimeError("Nessun evento valido dopo allineamento.")

mean = np.full(npts, np.nan, dtype=np.float64)
m = counts > 0
mean[m] = sum_bins[m] / counts[m]

# ---- salva TH1D
h_mean = ROOT.TH1D("h_mean_wvf", "Mean shifted normalized waveform;sample;mean(wvf/charge)", npts, 0, npts*bin_width)
for i in range(npts):
    if np.isfinite(mean[i]):
        h_mean.SetBinContent(i+1, mean[i])

outname = f"{plot_folder}/{label}_{run}_mean_wvf_shifted.root"
out = ROOT.TFile(outname, "RECREATE")
h_mean.Write()
h_base.Write()
h_charge_mid.Write()
h_charge_late.Write()
out.Close()
print("Salvato:", outname)










'''

#!/usr/bin/env python3
import os
import numpy as np
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)

base_folder = '/home/carla/Documents/ADC_data/Legend/'
label = '20260211_PDU_cold_test_44_7V'
run = 'run11'

plot_folder = '/var/www/html/202602_PDU_cold_test/'



f = ROOT.TFile.Open(os.path.join(base_folder, label, run, 'analyzed.root'))
if not f or f.IsZombie():
    raise RuntimeError(f"File non apribile: {path}")

t = f.Get('t1')
if not t:
    raise RuntimeError(f"Tree t1 non trovato")

# leggi solo quello che serve
t.SetBranchStatus("*", 0)
t.SetBranchStatus('Ch0.charge', 1)
t.SetBranchStatus('Ch0.wvf', 1)
t.SetBranchStatus('Ch0.fprompt', 1)

npts = 5000
q_buf = array('d', [0.0])
f_buf = array('d', [0.0])
w_buf = np.zeros(npts, dtype=np.float64)

# collega i buffer ai branch
t.SetBranchAddress('Ch0.charge', q_buf)
t.SetBranchAddress('Ch0.fprompt', f_buf)
t.SetBranchAddress('Ch0.wvf', w_buf)

sum_bins = np.zeros(npts, dtype=np.float64)
counts   = np.zeros(npts, dtype=np.int64)

n_tot = t.GetEntries()
n_used = 0

for i in range(n_tot):
    t.GetEntry(i)

    charge = float(q_buf[0])
    fprompt = float(f_buf[0])
    if not np.isfinite(charge) or charge == 0.0:
        continue
    if charge > 7e6:
        continue
    if fprompt > 0.4:
        continue
    
    w_norm = w_buf / charge
    valid = np.isfinite(w_norm)

    sum_bins[valid] += w_norm[valid]
    counts[valid]   += 1
    n_used += 1

print(f"Eventi totali: {n_tot}, usati: {n_used}")
if n_used == 0:
    raise RuntimeError("Nessun evento valido (charge=0 o non finito).")

mean = np.full(npts, np.nan, dtype=np.float64)
m = counts > 0
mean[m] = sum_bins[m] / counts[m]

h_mean = ROOT.TH1D("h_mean_wvf", "", npts, 0, npts)

for i in range(npts):
    if np.isfinite(mean[i]):
        h_mean.SetBinContent(i+1, mean[i])  # bin ROOT parte da 1

out = ROOT.TFile("%s/%s_%s_mean_wvf.root" %(plot_folder, label, run), "RECREATE")
h_mean.Write()
out.Close()
'''