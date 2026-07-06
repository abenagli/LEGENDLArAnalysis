#!/usr/bin/env python3
import numpy as np
import ROOT


def _th1_to_np(h, N, start_bin=1):
    a = np.empty(N, dtype=np.float64)
    for i in range(N):
        a[i] = h.GetBinContent(start_bin + i)
    return a


def _np_to_th1(a, name, title, x0, dx):
    """
    Crea TH1D con asse X in tempo.
    dx: bin width (qui in ns)
    x0: xmin (qui in ns)
    """
    N = len(a)
    hout = ROOT.TH1D(name, title, N, x0, x0 + N * dx)
    for i, v in enumerate(a, start=1):
        hout.SetBinContent(i, float(v))
    return hout


def _power_spectrum_rfft(x_pad):
    X = np.fft.rfft(x_pad)
    return (np.abs(X) ** 2)


def _np_to_th1_freq(power, name, title, dt_ns):
    """
    power: array len = L//2+1 (rfft)
    dt_ns: bin width in ns
    Output frequency in GHz (since 1/ns = GHz)
    """
    n = len(power)

    if dt_ns <= 0:
        raise ValueError(f"dt_ns non valido: {dt_ns}")

    fs_GHz = 1.0 / dt_ns       # GHz (samples per ns)
    fmax = fs_GHz / 2.0        # Nyquist
    df = fmax / (n - 1) if n > 1 else fmax

    # ultimo edge a fmax+df per includere il bin a f~fmax
    h = ROOT.TH1D(name, title, n, 0.0, fmax + df)
    for i, v in enumerate(power, start=1):
        h.SetBinContent(i, float(v))
    return h


def _next_pow2_at_least(n: int) -> int:
    L = 1
    while L < n:
        L *= 2
    return L


def _lowpass_rfft_brickwall(x: np.ndarray, dt_ns: float, fc_GHz: float, pad_to_pow2: bool = True) -> np.ndarray:
    """
    Low-pass "ideale" in frequenza (brick-wall):
      - rfft
      - azzera f > fc
      - irfft

    dt_ns: passo temporale in ns -> frequenze in GHz (1/ns)
    fc_GHz: frequenza di taglio (GHz)
    """
    if dt_ns <= 0:
        raise ValueError(f"dt_ns non valido: {dt_ns}")
    if fc_GHz <= 0:
        raise ValueError(f"fc_GHz non valido: {fc_GHz}")

    N = len(x)
    L = _next_pow2_at_least(N) if pad_to_pow2 else N

    x_pad = np.zeros(L, dtype=np.float64)
    x_pad[:N] = x

    X = np.fft.rfft(x_pad)
    freqs_GHz = np.fft.rfftfreq(L, d=dt_ns)  # 1/ns = GHz

    X[freqs_GHz > fc_GHz] = 0.0

    x_lp = np.fft.irfft(X, n=L)
    return np.real(x_lp[:N])


def deconvolve_remove_h1_from_h2_fft(
    h1, h2,
    N,
    eps=1e-12,
    pad_factor=2,
    start_bin=1,
    fc_lp_GHz=0.03,
    out_name="h_response",
    out_title="Deconvolution result;Time [ns];Amplitude",
):
    """
    Restituisce:
      - h_resp (TH1D)
      - h1_trunc, h2_trunc (TH1D)                 [raw]
      - h1_trunc_lp, h2_trunc_lp (TH1D)           [LP]
      - P_h1_raw, P_h2_raw (TH1D)                 [raw spectra]
      - P_h1_lp,  P_h2_lp  (TH1D)                 [LP spectra]
      - P_hresp (TH1D)                            [response spectrum]
    Deconvoluzione: H = Y * conj(X) / (|X|^2 + eps)
    """
    print(
        f"Deconvolution: remove {h1.GetName()} effect from {h2.GetName()}, "
        f"N={N}, eps={eps}, pad_factor={pad_factor}, start_bin={start_bin}, fc_lp={fc_lp_GHz} GHz"
    )

    # Assumo assi in ns
    dx1 = h1.GetXaxis().GetBinWidth(1)  # ns
    dx2 = h2.GetXaxis().GetBinWidth(1)  # ns
    if abs(dx1 - dx2) > 1e-12:
        raise ValueError(f"Bin width diverso: dx1={dx1} ns, dx2={dx2} ns")

    x0 = h1.GetXaxis().GetXmin() + (start_bin - 1) * dx1  # ns

    # segnali troncati (raw)
    x = _th1_to_np(h1, N, start_bin=start_bin)
    y = _th1_to_np(h2, N, start_bin=start_bin)

    h1_trunc = _np_to_th1(
        x, "h1_trunc",
        f"{h1.GetName()} truncated;Time [ns];Amplitude",
        x0=x0, dx=dx1
    )
    h2_trunc = _np_to_th1(
        y, "h2_trunc",
        f"{h2.GetName()} truncated;Time [ns];Amplitude",
        x0=x0, dx=dx1
    )

    # FFT length: almeno pad_factor*N e poi pow2
    target = int(pad_factor * N)
    L = _next_pow2_at_least(target)

    # ===== spettri RAW (prima del filtro) con stesso padding L =====
    x_raw_pad = np.zeros(L, dtype=np.float64)
    y_raw_pad = np.zeros(L, dtype=np.float64)
    x_raw_pad[:N] = x
    y_raw_pad[:N] = y

    P1_raw = _power_spectrum_rfft(x_raw_pad)
    P2_raw = _power_spectrum_rfft(y_raw_pad)

    hP1_raw = _np_to_th1_freq(
        P1_raw, "P_h1_raw",
        "Power spectrum |FFT(h1_raw)|^{2};Frequency [GHz];Power",
        dt_ns=dx1
    )
    hP2_raw = _np_to_th1_freq(
        P2_raw, "P_h2_raw",
        "Power spectrum |FFT(h2_raw)|^{2};Frequency [GHz];Power",
        dt_ns=dx1
    )

    # ===== LP FILTER prima della deconvoluzione =====
    x_lp = _lowpass_rfft_brickwall(x, dt_ns=dx1, fc_GHz=fc_lp_GHz, pad_to_pow2=True)
    y_lp = _lowpass_rfft_brickwall(y, dt_ns=dx1, fc_GHz=fc_lp_GHz, pad_to_pow2=True)

    h1_trunc_lp = _np_to_th1(
        x_lp, "h1_trunc_lp",
        f"{h1.GetName()} truncated (LP {fc_lp_GHz:.3f} GHz);Time [ns];Amplitude",
        x0=x0, dx=dx1
    )
    h2_trunc_lp = _np_to_th1(
        y_lp, "h2_trunc_lp",
        f"{h2.GetName()} truncated (LP {fc_lp_GHz:.3f} GHz);Time [ns];Amplitude",
        x0=x0, dx=dx1
    )

    # padding L coerente per LP + deconvoluzione
    x_pad = np.zeros(L, dtype=np.float64)
    y_pad = np.zeros(L, dtype=np.float64)
    x_pad[:N] = x_lp
    y_pad[:N] = y_lp

    # ===== spettri LP =====
    P1_lp = _power_spectrum_rfft(x_pad)
    P2_lp = _power_spectrum_rfft(y_pad)

    hP1_lp = _np_to_th1_freq(
        P1_lp, "P_h1_lp",
        "Power spectrum |FFT(h1_LP)|^{2};Frequency [GHz];Power",
        dt_ns=dx1
    )
    hP2_lp = _np_to_th1_freq(
        P2_lp, "P_h2_lp",
        "Power spectrum |FFT(h2_LP)|^{2};Frequency [GHz];Power",
        dt_ns=dx1
    )

    # ===== deconvoluzione su segnali LP =====
    X = np.fft.rfft(x_pad)
    Y = np.fft.rfft(y_pad)

    denom = (np.abs(X) ** 2 + eps)
    H = (Y * np.conj(X)) / denom

    h_time = np.fft.irfft(H, n=L)
    h = np.real(h_time[:N])

    h_resp = _np_to_th1(
        h, out_name, out_title,
        x0=x0, dx=dx1
    )

    # spettro della risposta
    h_pad = np.zeros(L, dtype=np.float64)
    h_pad[:N] = h
    Presp = _power_spectrum_rfft(h_pad)

    hPresp = _np_to_th1_freq(
        Presp, "P_hresp",
        "Power spectrum |FFT(h_response)|^{2};Frequency [GHz];Power",
        dt_ns=dx1
    )

    return (
        h_resp,
        h1_trunc, h2_trunc,
        h1_trunc_lp, h2_trunc_lp,
        hP1_raw, hP2_raw,
        hP1_lp,  hP2_lp,
        hPresp
    )


if __name__ == "__main__":
    # INPUT
    f1 = ROOT.TFile.Open("/var/www/html/202602_PDU_cold_test/20260211_PDU_cold_test_44_7V_run13_mean_wvf_shifted.root")
    f2 = ROOT.TFile.Open("/var/www/html/202602_PDU_cold_test/20260211_PDU_cold_test_44_7V_run11_mean_wvf_shifted.root")

    if not f1 or f1.IsZombie():
        raise RuntimeError("Errore apertura f1")
    if not f2 or f2.IsZombie():
        raise RuntimeError("Errore apertura f2")

    h1 = f1.Get("h_mean_wvf")   # effetto/response da rimuovere
    h2 = f2.Get("h_mean_wvf")   # segnale osservato
    if not h1 or not h2:
        raise RuntimeError("Non trovo 'h_mean_wvf' in uno dei file.")

    # Parametri
    N = 4096
    eps = 1e-4
    pad_factor = 1
    start_bin = 1
    fc_lp_GHz = 0.03

    (
        hresp,
        h1_trunc, h2_trunc,
        h1_trunc_lp, h2_trunc_lp,
        P1_raw, P2_raw,
        P1_lp,  P2_lp,
        Presp
    ) = deconvolve_remove_h1_from_h2_fft(
        h1, h2,
        N=N,
        eps=eps,
        pad_factor=pad_factor,
        start_bin=start_bin,
        fc_lp_GHz=fc_lp_GHz,
        out_name="h_response",
        out_title="Deconvolution result;Time [ns];Amplitude",
    )

    # OUTPUT
    out_path = "/var/www/html/202602_PDU_cold_test/20260211_PDU_cold_test_44_7V_decon.root"
    out = ROOT.TFile(out_path, "RECREATE")

    # time domain
    h1_trunc.Write()
    h2_trunc.Write()
    h1_trunc_lp.Write()
    h2_trunc_lp.Write()
    hresp.Write()

    # freq domain (RAW + LP + response)
    P1_raw.Write()
    P2_raw.Write()
    P1_lp.Write()
    P2_lp.Write()
    Presp.Write()

    out.Close()

    print("Scritto:", out_path)
    print("Oggetti salvati:")
    print("  h1_trunc, h2_trunc, h1_trunc_lp, h2_trunc_lp, h_response")
    print("  P_h1_raw, P_h2_raw, P_h1_lp, P_h2_lp, P_hresp")