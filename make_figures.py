# SPDX-License-Identifier: MIT OR Apache-2.0
# Copyright (c) 2026 The OpenKeito contributors
"""Make the IEEE 39-bus figures and tables of the paper (F14-F18, T8-T10, OpenKeito results) and their numbers.

    python make_figures.py [--runs runs] [--out figures]

Inputs: the CIM (``cim/``), the models (``models/``) and the OpenKeito results of ``run_cases.py``
(``runs/``).  Outputs in ``figures/``:

* ``f14_oneline.pdf``, ``f15_steady_state.{pdf,csv}``, ``f16_pv_startup.{pdf,csv}``,
  ``f17_fault_current.{pdf,csv}``, ``f17_fault_current_cycles.csv``, ``f18_response_20s.{pdf,csv}``;
* ``t08_mapping.tex``, ``t09_sequence_impedance.tex``, ``t10_sensitivity.tex`` (booktabs tabular bodies);
* ``numbers.json`` (every number quoted in the text) and ``manifest.csv`` (OpenKeito version, sha256 of the
  CIM inputs and of every output).

The same inputs give the same files: PDFs carry no dates, CSV files are UTF-8 with LF and ten significant
digits.  ``numbers.json`` has ``openkeito_checks.compute_time`` only when ``runs/timing.json``
(``run_cases.py --timing``) exists; computing time depends on the machine.
"""
import argparse
import cmath
import csv
import hashlib
import json
import math
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cim_to_openkeito as C  # noqa: E402

CIM, IC = HERE / 'cim' / 'IEEE39.xml', HERE / 'cim' / 'IEEE39_ic.xml'
MODELS = HERE / 'models'
OPENKEITO = 'https://github.com/yutanakamurajp/OpenKeito'
W_SINGLE, W_DOUBLE = 3.358, 6.989          # IEEE Access column widths (in)
COLORS = ['#0072B2', '#D55E00', '#009E73', '#AA4499']
LINES = ['-', '--', '-.', ':']
MARKERS = ['o', 's', '^', 'D']
F, WW = 60.0, 2 * math.pi * 60.0
SR = 480.0                                  # PV plant rating (MVA)
FAULT2 = C.CASES['slg16']['fault']
# F14: the buses lie at the CIM diagram coordinates (TextDiagramObject).  Buses that overlap there (each
# generator bus sits on the other bus of its transformer; buses 10-13, 19/20/33/34, 21-24) are moved by
# these amounts (diagram units, x right, y up), and bus numbers are placed by these offsets (points).
F14_BUS_SHIFT = {'30': (-0.053, -0.051), '37': (-0.051, -0.046), '38': (-0.047, -0.047), '31': (0.093, 0.1),
                 '6': (0.036, 0.052), '10': (-0.023, 0.067), '11': (0.053, 0.048), '12': (-0.003, -0.014),
                 '13': (-0.06, 0.044), '32': (-0.024, 0.134), '21': (-0.001, 0.047), '24': (-0.024, -0.051),
                 '35': (-0.042, 0.056), '36': (-0.049, -0.049), '19': (0.024, 0.026), '20': (-0.023, 0.024),
                 '33': (-0.029, -0.021), '34': (0.028, -0.029)}
F14_LABEL_OFFSET_PT = {'default': (2.5, 2.5), '30': (-11, -7), '37': (-11, -7), '38': (-11, -7), '33': (-11, -7),
                       '34': (3, -7), '20': (-9, 3), '19': (3, 2), '6': (3, -6), '14': (-3, -9), '5': (-3, -8),
                       '10': (-9, 0), '11': (3, -3), '12': (-3, -8), '13': (-9, 1), '32': (3, 0), '23': (-9, 2),
                       '35': (-10, 0), '36': (-11, -4), '24': (-3, -8)}
RUNS = HERE / 'runs'


def style():
    plt.rcParams.update({
        'font.family': ['DejaVu Sans'], 'mathtext.fontset': 'dejavusans', 'font.size': 8.0, 'axes.labelsize': 8.0,
        'axes.titlesize': 8.0, 'xtick.labelsize': 8.0, 'ytick.labelsize': 8.0, 'legend.fontsize': 7.0,
        'legend.frameon': False, 'legend.handlelength': 2.6, 'lines.linewidth': 1.2, 'lines.markersize': 4.5,
        'axes.linewidth': 0.6, 'axes.grid': True, 'axes.axisbelow': True, 'grid.color': '#d9d9d9', 'grid.linewidth': 0.4,
        'xtick.direction': 'in', 'ytick.direction': 'in', 'figure.dpi': 100, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'figure.constrained_layout.use': True, 'svg.hashsalt': 'ieee39'})


def ser(k, marker=False, line=True):
    s = {'color': COLORS[k], 'linestyle': LINES[k] if line else 'none'}
    if marker:
        s.update(marker=MARKERS[k], markerfacecolor='white', markeredgecolor=COLORS[k], markeredgewidth=1.0)
    return s


def save(fig, path):
    fig.savefig(path, metadata={'CreationDate': None, 'ModDate': None, 'Producer': None, 'Creator': None})
    plt.close(fig)


def wcsv(path, df):
    def f(v):
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return ''
        if isinstance(v, (float, np.floating)):
            return repr(float(f'{float(v):.10g}'))
        return str(v)
    with open(path, 'w', newline='', encoding='utf-8') as fp:
        w = csv.writer(fp, lineterminator='\n')
        w.writerow(df.columns)
        for r in df.itertuples(index=False):
            w.writerow([f(v) for v in r])


# ---------------------------------------------------------------- results

def standard(case):
    """runs/<case>/ok.csv: the recorded signals of result.csv under short names and in MW, Mvar, kA, pu."""
    d = RUNS / case
    r = pd.read_csv(d / 'run' / 'result.csv')
    o = pd.DataFrame({'t_s': r['time_s']})
    for s in C.RECORD_PV:
        o[f'IBR_{s}'] = r[f'S(IBR30_1_{s})']
    o['IBR_P_pu'], o['IBR_Q_pu'] = o['IBR_PDQ0'], o['IBR_QDQ0']
    o['IBR_P_MW'], o['IBR_Q_Mvar'] = o['IBR_PDQ0'] * SR, o['IBR_QDQ0'] * SR
    o['IBR_F_pu'] = o['IBR_WE'] / C.WC
    f = C.CASES[case].get('fault')
    if f:
        o['IF_kA'] = r[f"I(F{f['bus']}.{f['phases'][-1]})"] / 1e3
        for p in 'abc':
            o[f"V{f['bus']}_{p.upper()}_V"] = r[f"V(B{f['bus']}.{p})"]
    o['SW_P_MW'] = r['S(SWP_W)']
    for k in range(32, 40):
        g = f'G{k}_1'
        o[f'G{k}_F_pu'], o[f'G{k}_V_pu'] = r[f'S({g}_OMG)'], r[f'S({g}_EA)']
        o[f'G{k}_Pm_pu'], o[f'G{k}_Efd_pu'], o[f'G{k}_P_pu'] = r[f'S({g}_TMO)'], r[f'S({g}_EFO)'], r[f'S({g}_P)']
    o.to_csv(d / 'ok.csv', index=False, encoding='utf-8')


def rd(case):
    return pd.read_csv(RUNS / case / 'ok.csv', encoding='utf-8')


def mean(d, col, w):
    return float(d[(d.t_s >= w[0]) & (d.t_s < w[1])][col].mean())


def peak(d, col, w, fn):
    s = d[(d.t_s >= w[0]) & (d.t_s <= w[1])]
    i = s[col].idxmax() if fn == 'max' else s[col].idxmin()
    return float(s.at[i, 't_s']), float(s.at[i, col])


def cycles(t, y, t0, n):
    """Half the peak-to-peak (AC) and the mid value (DC) of each cycle after t0."""
    out = []
    for k in range(n):
        m = (t >= t0 + k / F - 1e-9) & (t < t0 + (k + 1) / F - 1e-9)
        out.append(((y[m].max() - y[m].min()) / 2, (y[m].max() + y[m].min()) / 2))
    return out


def phase(t, y, t0):
    """Phase angle (cosine reference) at t0 of the 60 Hz fit over the cycle before t0."""
    m = (t >= t0 - 1 / F) & (t < t0)
    X = np.column_stack([np.cos(WW * t[m]), np.sin(WW * t[m]), np.ones(m.sum())])
    c = np.linalg.lstsq(X, y[m], rcond=None)[0]
    return math.degrees(cmath.phase(complex(c[0], -c[1]) * np.exp(1j * WW * t0)))


def thin(d, step_ms=1.0):
    ms = np.round(d.t_s.to_numpy() * 1000, 6)
    return d[np.abs(ms / step_ms - np.round(ms / step_ms)) < 1e-6].reset_index(drop=True)


def net():
    return C.network(C.load_cim(CIM, IC))


def steady(n):
    """Phasor solution (nofault) against the CIM initial conditions: buses and machines."""
    p = pd.read_csv(RUNS / 'nofault' / 'run' / 'phasors.csv')
    ph = {q: complex(a, b) for q, a, b in zip(p.quantity, p.real, p.imag)}
    rows = []
    for name in sorted(n['buses'], key=int):
        v0, a0 = n['ic'][name]
        v = ph[f'V(B{name}.a)']
        vm = abs(v) / (n['buses'][name]['kv'] * 1e3 / math.sqrt(3))
        va = math.degrees(cmath.phase(v)) - n['shift'][name]
        rows.append({'element': f'bus {name}', 'quantity': 'v_pu', 'cim_ic': v0, 'openkeito': vm, 'diff': vm - v0})
        rows.append({'element': f'bus {name}', 'quantity': 'angle_deg', 'cim_ic': a0, 'openkeito': va,
                     'diff': (va - a0 + 180) % 360 - 180})
    for g in sorted(n['gens'] + n['pvs'], key=lambda g: int(g['bus'])):
        cid = f"SW{g['name']}" if g['bus'] == C.SWING_BUS else (f"IBRS{g['name']}" if g in n['pvs'] else f"G{g['name']}")
        s = -3 * ph[f"V(B{g['bus']}.a)"] * ph[f'I({cid}.a)'].conjugate() / 1e6
        rows.append({'element': f"gen {g['bus']}", 'quantity': 'p_mw', 'cim_ic': g['p_mw'], 'openkeito': s.real,
                     'diff': s.real - g['p_mw']})
        rows.append({'element': f"gen {g['bus']}", 'quantity': 'q_mvar', 'cim_ic': g['q_mvar'], 'openkeito': s.imag,
                     'diff': s.imag - g['q_mvar']})
    return pd.DataFrame(rows)


def start_deviation(case, t_end=0.01):
    """Largest deviation before 0.01 s (when the WECC P and Q references start to rise) of the recorded voltages and
    currents from the initial phasor solution, relative to the peak value."""
    p = pd.read_csv(RUNS / case / 'run' / 'phasors.csv')
    ph = {q: complex(a, b) for q, a, b in zip(p.quantity, p.real, p.imag)}
    r = pd.read_csv(RUNS / case / 'run' / 'result.csv')
    m = r.time_s < t_end
    t = r.time_s[m].to_numpy()
    worst = 0.0
    for c in r.columns:
        if c in ph:
            x = math.sqrt(2) * np.real(ph[c] * np.exp(1j * WW * t))
            worst = max(worst, float(np.abs(r[c][m].to_numpy() - x).max() / (math.sqrt(2) * abs(ph[c]) + 1e-9)))
    return worst


def oscillation(t, y, t0, win=0.005, skip=2e-4):
    """Numerical-oscillation index: largest second difference over 5 ms from 0.2 ms after an event, divided by the
    largest value within one cycle around the event."""
    m = (t >= t0 + skip) & (t < t0 + skip + win)
    yy = y[m]
    d2 = np.abs(yy[1:-1] - (yy[:-2] + yy[2:]) / 2)
    return float(d2.max() / max(np.abs(y[(t >= t0 - 1 / F) & (t < t0 + 1 / F)]).max(), 1e-9))


def extra_numbers():
    d = rd('slg16')
    t = d.t_s.to_numpy()
    osc = {}
    for ev, t0 in (('inception', FAULT2['t_on']), ('clearing', FAULT2['t_off'])):
        for col, key in (('IF_kA', 'fault_current'), ('V16_A_V', 'v16_a'), ('V16_C_V', 'v16_c')):
            osc[f'{ev}_{key}'] = oscillation(t, d[col].to_numpy(), t0)
    slip = {}
    for case in ('slg16', 'tpf02', 'slg16_20s'):
        x = rd(case)
        slip[case] = float(100 * max((x[f'G{k}_F_pu'] - 1).abs().max() for k in range(32, 40)))
    out = {'start_deviation_to_0p01s_rel': {c: start_deviation(c) for c in ('nofault', 'nofault_qreg0', 'slg16', 'tpf02')},
           'oscillation_index_slg16': osc, 'machine_slip_max_pct': slip,
           'definitions': {'start_deviation_to_0p01s_rel': 'max over recorded voltages and currents of max|y(t) - sqrt(2) Re(Y e^{j omega t})| for t < 0.01 s, divided by the peak sqrt(2)|Y|, where Y is the initial phasor solution (phasors.csv)',
                           'oscillation_index_slg16': 'max |y[k] - (y[k-1] + y[k+1])/2| over 5 ms from 0.2 ms after the event, divided by max|y| within one cycle around the event',
                           'machine_slip_max_pct': 'max over machines G32-G39 and time of |speed - 1| in percent',
                           'compute_time': 'OpenKeito alone (openkeito run, wall time in summary.json), model building and file transfer excluded'}}
    tp = RUNS / 'timing.json'
    if tp.exists():
        tm = json.loads(tp.read_text(encoding='utf-8'))
        out['compute_time'] = {'machine': tm['cpu'], 'os': tm['os'], 'condition': f"{tm['affinity']}, threads {tm['threads']}, "
                               f"{tm['statistic']} ({tm['repeat']} runs)",
                               'cases': {k: {'elapsed_s': v['elapsed_s'], 'us_per_step': v['us_per_step'], 'steps': v['steps']}
                                         for k, v in tm['cases'].items()}}
    return out


# ---------------------------------------------------------------- figures

def f14(n, out):
    """One-line diagram from the CIM diagram coordinates (TextDiagramObject of every bus)."""
    ns = {'cim': 'http://www.ucaiug.org/grid18v15#', 'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'}
    root = ET.parse(CIM).getroot()
    rdf = '{%s}' % ns['rdf']
    obj = {}
    for o in root.findall('cim:TextDiagramObject', ns):
        obj[o.get(rdf + 'about')] = o.find('cim:IdentifiedObject.name', ns).text
    pos = {}
    for pnt in root.findall('cim:DiagramObjectPoint', ns):
        oid = pnt.find('cim:DiagramObjectPoint.DiagramObject', ns).get(rdf + 'resource')
        if oid in obj:
            pos[obj[oid]] = (float(pnt.find('cim:DiagramObjectPoint.xPosition', ns).text),
                             float(pnt.find('cim:DiagramObjectPoint.yPosition', ns).text))
    for b, (dx, dy) in F14_BUS_SHIFT.items():
        pos[b] = (pos[b][0] + dx, pos[b][1] + dy)
    fig, ax = plt.subplots(figsize=(W_SINGLE, 3.4))
    ax.grid(False)
    for ln in n['lines']:
        (x1, y1), (x2, y2) = pos[ln['from']], pos[ln['to']]
        ax.plot([x1, x2], [y1, y2], color='0.35', lw=0.7)
    for t in n['xfmrs']:
        (x1, y1), (x2, y2) = pos[t['from']], pos[t['to']]
        ax.plot([x1, x2], [y1, y2], color='0.35', lw=0.7, ls=(0, (2, 1)))
    gen_bus = {g['bus'] for g in n['gens']}
    for b, (x, y) in sorted(pos.items(), key=lambda kv: int(kv[0])):
        if b == '30':
            ax.plot(x, y, marker='s', ms=6, color=COLORS[2], mec='k', mew=0.5, zorder=3)
        elif b in gen_bus:
            ax.plot(x, y, marker='o', ms=5, color='white', mec='k', mew=0.7, zorder=3)
        else:
            ax.plot(x, y, marker='o', ms=2.5, color='k', zorder=3)
        if b in (FAULT2['bus'], '2'):
            ax.plot(x, y, marker='x', ms=8, color=COLORS[1], mew=1.5, zorder=4)
        ax.annotate(b, (x, y), xytext=F14_LABEL_OFFSET_PT.get(b, F14_LABEL_OFFSET_PT['default']),
                    textcoords='offset points', fontsize=5.5)
    h = [plt.Line2D([], [], marker='s', ms=6, color=COLORS[2], mec='k', ls='none'),
         plt.Line2D([], [], marker='o', ms=5, color='white', mec='k', ls='none'),
         plt.Line2D([], [], marker='x', ms=7, color=COLORS[1], mew=1.5, ls='none'),
         plt.Line2D([], [], color='0.35', lw=0.7), plt.Line2D([], [], color='0.35', lw=0.7, ls=(0, (2, 1)))]
    ax.legend(h, ['PV plant', 'Synchronous machine', 'Fault', 'Line', 'Transformer'], loc='upper center',
              bbox_to_anchor=(0.5, -0.01), ncols=3, fontsize=6)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal', adjustable='datalim')
    for s in ax.spines.values():
        s.set_visible(False)
    save(fig, out / 'f14_oneline.pdf')


def f15(st, out):
    fig, axs = plt.subplots(2, 1, figsize=(W_SINGLE, 3.4))
    v = st[st.quantity == 'v_pu']
    buses = [int(e.split()[1]) for e in v.element]
    axs[0].semilogy(buses, v['diff'].abs().clip(lower=1e-9), **ser(0, marker=True, line=False), label='OpenKeito')
    g = st[st.quantity == 'p_mw']
    gb = [int(e.split()[1]) for e in g.element]
    axs[1].semilogy(gb, g['diff'].abs().clip(lower=1e-9), **ser(0, marker=True, line=False), label='OpenKeito')
    axs[0].set_xlabel('Bus'); axs[0].set_ylabel('|ΔV| (pu)')
    axs[1].set_xlabel('Generator bus'); axs[1].set_ylabel('|ΔP| (MW)'); axs[1].set_xticks(gb)
    save(fig, out / 'f15_steady_state.pdf')


def f16(out):
    d = thin(rd('nofault'))
    d = d[d.t_s <= 5.0]
    fig, axs = plt.subplots(3, 1, figsize=(W_SINGLE, 3.8), sharex=True)
    axs[0].plot(d.t_s, d.IBR_P_MW, **ser(0), label='OpenKeito')
    axs[1].plot(d.t_s, d.IBR_Q_Mvar, **ser(0))
    axs[2].plot(d.t_s, d.IBR_F_pu, **ser(0))
    axs[0].set_ylabel('P (MW)'); axs[1].set_ylabel('Q (Mvar)'); axs[2].set_ylabel('PLL freq. (pu)')
    axs[0].set_ylim(0, 300); axs[1].set_ylim(-5, 110); axs[2].set_ylim(0.995, 1.005)
    axs[2].set_xlabel('Time (s)'); axs[2].set_xlim(0, 5)
    save(fig, out / 'f16_pv_startup.pdf')
    return pd.DataFrame({'time_s': d.t_s, 'p_mw': d.IBR_P_MW, 'q_mvar': d.IBR_Q_Mvar, 'pll_freq_pu': d.IBR_F_pu})


def f17(out):
    d = rd('slg16')
    t0, t1 = FAULT2['t_on'], FAULT2['t_off']
    n = 12
    cy = cycles(d.t_s.to_numpy(), d.IF_kA.to_numpy(), t0, n)
    fig, axs = plt.subplots(2, 1, figsize=(W_SINGLE, 3.4))
    w = d[(d.t_s >= t0 - 0.02) & (d.t_s <= t1 + 0.05)]
    axs[0].plot(w.t_s, w.IF_kA, **ser(0), lw=0.9, label='OpenKeito')
    axs[0].set_xlabel('Time (s)'); axs[0].set_ylabel('Fault current (kA)'); axs[0].set_xlim(t0 - 0.02, t1 + 0.05)
    k = np.arange(1, n + 1)
    axs[1].plot(k, [c[0] for c in cy], **ser(0, marker=True), lw=0.8, label='AC, OpenKeito')
    axs[1].plot(k, [c[1] for c in cy], **dict(ser(0, marker=True), linestyle='--'), lw=0.8, label='DC, OpenKeito')
    axs[1].axhline(0, color='0.6', lw=0.4)
    axs[1].set_xlabel('Cycle after inception'); axs[1].set_ylabel('AC / DC (kA)'); axs[1].set_xticks(k)
    axs[1].legend(loc='center right', fontsize=6)
    save(fig, out / 'f17_fault_current.pdf')
    wf = thin(w, 0.1)
    wave = pd.DataFrame({'time_s': wf.t_s, 'fault_current_ka': wf.IF_kA})
    cyc = pd.DataFrame({'cycle': k, 'ac_ka': [c[0] for c in cy], 'dc_ka': [c[1] for c in cy]})
    return wave, cyc, phase(d.t_s.to_numpy(), d.V16_C_V.to_numpy(), t0)


def f18(out):
    d = rd('slg16_20s')
    d10 = d[(np.abs(np.round(d.t_s * 1000) % 10) < 1e-6) | ((d.t_s >= 11.9) & (d.t_s <= 12.6))]
    fig, axs = plt.subplots(2, 1, figsize=(W_DOUBLE, 3.6))
    axs[0].plot(d10.t_s, d10.IBR_Q_Mvar, **ser(0), lw=0.9, label='OpenKeito')
    z = d10[(d10.t_s >= 11.8) & (d10.t_s <= 16)]
    axs[1].plot(z.t_s, z.G35_F_pu, **ser(0), lw=0.9, label='OpenKeito')
    axs[0].set_xlim(0, 20); axs[0].set_ylim(-10, 180); axs[0].set_xlabel('Time (s)'); axs[0].set_ylabel('PV Q (Mvar)')
    axs[1].set_xlim(11.8, 16); axs[1].set_ylim(0.994, 1.009); axs[1].set_xlabel('Time (s)'); axs[1].set_ylabel('G35 speed (pu)')
    save(fig, out / 'f18_response_20s.pdf')
    return pd.DataFrame({'time_s': d10.t_s, 'pv_q_mvar': d10.IBR_Q_Mvar, 'g35_speed_pu': d10.G35_F_pu})


# ---------------------------------------------------------------- tables

def tex(path, head, rows, align, note=None):
    lines = ['% Made by make_figures.py'] + ([f'% {note}'] if note else []) + [f'\\begin{{tabular}}{{{align}}}', '\\toprule', ' & '.join(head) + r' \\', '\\midrule']
    lines += [' & '.join(str(c) for c in r) + r' \\' for r in rows]
    lines += ['\\bottomrule', '\\end{tabular}', '']
    path.write_text('\n'.join(lines), encoding='utf-8', newline='\n')


def t08(out):
    m = json.loads((MODELS / 'slg16_20s.json').read_text(encoding='utf-8'))
    cnt = {}
    for c in m['components']:
        key = c['type'] + (' (magnetizing)' if c['name'].startswith('MAG') else '')
        cnt[key] = cnt.get(key, 0) + 1
    rows = [('ACLineSegment', r'\texttt{bergeron\_line}', cnt['bergeron_line']),
            ('PowerTransformer', r'\texttt{three\_phase\_transformer}', cnt['three_phase_transformer']),
            ('TransformerCoreAdmittance', r'\texttt{three\_phase\_load} (excluded from the phasor solution)',
             cnt['three_phase_load (magnetizing)']),
            ('EnergyConsumer', r'\texttt{three\_phase\_load}', cnt['three_phase_load']),
            ('SynchronousMachine (bus 32--39)', r'\texttt{synchronous\_machine} + exciter and governor blocks',
             cnt['synchronous_machine']),
            ('SynchronousMachine (bus 31, swing)', r'\texttt{three\_phase\_source} behind $R_a + jX_d^\prime$', 1),
            ('PowerElectronicsConnection (PV)', r'\texttt{three\_phase\_controlled\_current\_source} + REGC\_A/REEC\_A/REPC\_A', 1),
            ('(initializing source of the PV)', r'\texttt{three\_phase\_source} + 3 \texttt{switch}', 1),
            ('(fault)', r'\texttt{fault}', cnt['fault'])]
    nf = json.loads((MODELS / 'nofault.json').read_text(encoding='utf-8'))
    rows.append(('Total (no-fault case)', f"{len(nf['components'])} components, {len(nf['controls'])} control blocks", ''))
    tex(out / 't08_mapping.tex', ['CIM class', 'OpenKeito', 'Count'], rows, 'llr')
    return {'components_nf': len(nf['components']), 'controls_nf': len(nf['controls']),
            'components_20s': len(m['components'])}


def scan_table():
    """runs/scan.csv: the sequence impedances of openkeito scan (runs/scan/<name>/scan.csv) against the CIM."""
    exp = json.loads((MODELS / 'scan' / 'expected_ohm.json').read_text(encoding='utf-8'))
    rows = []
    for n, e in exp.items():
        sc = pd.read_csv(RUNS / 'scan' / n / 'scan.csv')
        for k, (re_, im) in e.items():
            zs = cmath.rect(float(sc[f'{k}mag(P)'][0]), math.radians(float(sc[f'{k}deg(P)'][0])))
            ze = complex(re_, im)
            rows.append({'test': n, 'seq': k, 'r_ohm': zs.real, 'x_ohm': zs.imag, 'r_cim_ohm': ze.real, 'x_cim_ohm': ze.imag,
                         'z_rel_diff': abs(abs(zs) / abs(ze) - 1), 'x_rel_diff': abs(zs.imag / ze.imag - 1)})
    pd.DataFrame(rows).to_csv(RUNS / 'scan.csv', index=False, encoding='utf-8')
    return pd.read_csv(RUNS / 'scan.csv', encoding='utf-8')     # the table is read back from the file


def t09(out):
    s = scan_table()
    names = {'line_15_16_1_short': 'Line 15--16, far end shorted', 'line_15_16_1_open': 'Line 15--16, far end open',
             'line_16_19_1_short': 'Line 16--19, far end shorted', 'line_16_19_1_open': 'Line 16--19, far end open',
             'xfmr_19_33_0_1_short': 'Transformer 19--33 (YNd1), LV shorted', 'xfmr_19_33_0_1_open': 'Transformer 19--33, LV open',
             'xfmr_19_20_0_1_short': 'Transformer 19--20 (YNyn0), LV shorted'}
    rows = []
    for _, r in s.iterrows():
        if not (r.test in names or r.test in ('gen_32_1', 'gen_39_1', 'swing_31_1')):
            continue
        lab = names.get(r.test, {'gen_32_1': 'Machine G32', 'gen_39_1': 'Machine G39', 'swing_31_1': 'Swing source (bus 31)'}.get(r.test))
        rows.append((lab, {'Z1': 'pos.', 'Z2': 'neg.', 'Z0': 'zero'}[r.seq], f'{r.r_ohm:.4g} + j{r.x_ohm:.5g}',
                     f'{r.r_cim_ohm:.4g} + j{r.x_cim_ohm:.5g}', f'{r.z_rel_diff:.1e}'))
    rows.append((f'All {len(s)} items (2 lines, 2 transformers, 9 machines, swing source)', '', '', '',
                 f'max {s.z_rel_diff.max():.1e}'))
    tex(out / 't09_sequence_impedance.tex', ['Element', 'Seq.', r'Scan ($\Omega$)', r'CIM ($\Omega$)', r'$|\Delta Z|/|Z|$'],
        rows, 'llrrr')
    return {'items': int(len(s)), 'max_rel_diff_z': float(s.z_rel_diff.max()), 'max_rel_diff_x': float(s.x_rel_diff.max())}


def sens_values():
    base, okpi = rd('nofault'), rd('nofault_conditional_pi')
    out = {}
    for k, d in (('state', base), ('conditional', okpi)):
        t, q = peak(d, 'IBR_Q_Mvar', (0.3, 5.0), 'max')
        out[f'pi_{k}'] = {'q_peak_mvar': q, 'q_peak_time_s': t, 'q_4p5_5_mvar': mean(d, 'IBR_Q_Mvar', (4.5, 5.0))}
    for k, case in (('exclude', 'slg16_20s'), ('shunt_in_phasor', 'slg16_20s_mag_shunt'), ('transformer', 'slg16_20s_mag_transformer')):
        d = rd(case)
        t, q = peak(d, 'IBR_Q_Mvar', (0.3, 11.9), 'max')
        out[f'mag_{k}'] = {'q_prefault_mvar': mean(d, 'IBR_Q_Mvar', (11.5, 12.0)), 'q_peak_mvar': q, 'q_peak_time_s': t,
                           'q_20s_mvar': float(d.IBR_Q_Mvar.iloc[-20:].mean())}
    return out


def t10(out, sv):
    p, m = sv, sv
    rows = [('PI anti-windup (REEC\\_A, REPC\\_A)', 'state (integral held)', f"{p['pi_state']['q_peak_mvar']:.1f} at {p['pi_state']['q_peak_time_s']:.3f} s",
             f"{p['pi_state']['q_4p5_5_mvar']:.2f}"),
            ('', 'conditional integration', f"{p['pi_conditional']['q_peak_mvar']:.1f} at {p['pi_conditional']['q_peak_time_s']:.3f} s",
             f"{p['pi_conditional']['q_4p5_5_mvar']:.2f}"),
            ('\\midrule Magnetizing branches', 'excluded, enter at $t=0^+$', f"{m['mag_exclude']['q_peak_mvar']:.1f} at {m['mag_exclude']['q_peak_time_s']:.2f} s",
             f"{m['mag_exclude']['q_prefault_mvar']:.2f}"),
            ('', 'shunt in the phasor solution', f"{m['mag_shunt_in_phasor']['q_peak_mvar']:.1f} at {m['mag_shunt_in_phasor']['q_peak_time_s']:.2f} s",
             f"{m['mag_shunt_in_phasor']['q_prefault_mvar']:.2f}"),
            ('', 'inside the transformer', f"{m['mag_transformer']['q_peak_mvar']:.1f} at {m['mag_transformer']['q_peak_time_s']:.2f} s",
             f"{m['mag_transformer']['q_prefault_mvar']:.2f}")]
    tex(out / 't10_sensitivity.tex', ['Choice', 'Variant', 'PV Q peak after start-up (Mvar)', 'PV Q (Mvar)$^*$'], rows, 'llrr',
        note='* PI rows: mean over 4.5-5 s of the no-fault case; magnetizing rows: mean over 11.5-12.0 s of the 20 s case')


# ---------------------------------------------------------------- main

def manifest(out, version):
    rows = [('openkeito', f"{version} ({OPENKEITO})", '')]
    for p in (CIM, IC):
        rows.append((p.name, 'input (emthub CIM, Apache-2.0)', hashlib.sha256(p.read_bytes()).hexdigest()))
    for p in sorted(out.iterdir()):
        if p.name != 'manifest.csv':
            rows.append((p.name, 'output', hashlib.sha256(p.read_bytes()).hexdigest()))
    with open(out / 'manifest.csv', 'w', newline='', encoding='utf-8') as fp:
        w = csv.writer(fp, lineterminator='\n')
        w.writerow(['file', 'role', 'sha256'])
        w.writerows(rows)


def build(out):
    out.mkdir(parents=True, exist_ok=True)
    style()
    for case in C.CASES:
        standard(case)
    n = net()
    st = steady(n)
    nums = {}
    f14(n, out)
    f15(st, out)
    wcsv(out / 'f15_steady_state.csv', st)
    v, g = st[st.quantity == 'v_pu'], st[st.quantity == 'p_mw']
    nums['f15'] = {'max_abs_dv_pu': float(v['diff'].abs().max()),
                   'max_abs_dangle_deg': float(st[st.quantity == 'angle_deg']['diff'].abs().max()),
                   'max_abs_dp_mw': float(g['diff'].abs().max()),
                   'max_abs_dq_mvar': float(st[st.quantity == 'q_mvar']['diff'].abs().max())}
    wcsv(out / 'f16_pv_startup.csv', f16(out))
    nf = rd('nofault')
    t, q = peak(nf, 'IBR_Q_Mvar', (0.3, 5.0), 'max')
    we = nf[(nf.t_s >= 5.0) & (nf.t_s <= 6.0)].IBR_F_pu
    nums['f16'] = {'p_4p5_5_mw': mean(nf, 'IBR_P_MW', (4.5, 5.0)), 'q_4p5_5_mvar': mean(nf, 'IBR_Q_Mvar', (4.5, 5.0)),
                   'q_peak_mvar': q, 'q_peak_time_s': t, 'pll_freq_ptp_5_6_pu': float(we.max() - we.min())}
    na = rd('nofault_qreg0')
    nums['f16']['qreg0_zero'] = {'p_4p5_5_mw': mean(na, 'IBR_P_MW', (4.5, 5.0)), 'q_4p5_5_mvar': mean(na, 'IBR_Q_Mvar', (4.5, 5.0))}
    wave, cyc, ang = f17(out)
    wcsv(out / 'f17_fault_current.csv', wave)
    wcsv(out / 'f17_fault_current_cycles.csv', cyc)
    nums['f17'] = {'inception_angle_deg': ang, 'ac1_ka': float(cyc.ac_ka[0]), 'dc1_ka': float(cyc.dc_ka[0]),
                   'ac12_ka': float(cyc.ac_ka.iloc[-1]), 'dc12_ka': float(cyc.dc_ka.iloc[-1])}
    wcsv(out / 'f18_response_20s.csv', f18(out))
    d20 = rd('slg16_20s')
    nums['f18'] = {'q_prefault_mvar': mean(d20, 'IBR_Q_Mvar', (11.5, 12.0)),
                   'p_prefault_mw': mean(d20, 'IBR_P_MW', (11.5, 12.0)),
                   'q_peak_start': peak(d20, 'IBR_Q_Mvar', (0.3, 11.9), 'max'),
                   'q_peak_after_clearing': peak(d20, 'IBR_Q_Mvar', (12.2, 13.5), 'max'),
                   'q_valley': peak(d20, 'IBR_Q_Mvar', (13.0, 16.0), 'min'),
                   'q_20s_mvar': float(d20.loc[(d20.t_s - 20).abs().idxmin(), 'IBR_Q_Mvar']),
                   'fault_current_min_max_ka': [float(d20[(d20.t_s >= 12.0) & (d20.t_s <= 12.2)].IF_kA.min()),
                                                float(d20[(d20.t_s >= 12.0) & (d20.t_s <= 12.2)].IF_kA.max())],
                   'g35_speed_max': peak(d20, 'G35_F_pu', (12.0, 13.0), 'max'),
                   'g35_speed_min': peak(d20, 'G35_F_pu', (12.0, 13.0), 'min')}
    f3 = rd('tpf02')
    m = (f3.t_s >= 3.0) & (f3.t_s <= 3.1)
    nums['f02_3lg'] = {'vgpu_min': float(f3[m].IBR_VGPU.min()), 'vsup_max': float(f3[m].IBR_VSUP.max()),
                       'iqcmd_max': float(f3[m].IBR_IQCMD.max()), 'p_min_mw': float(f3[m].IBR_P_MW.min())}
    nums['t08'] = t08(out)
    nums['t09'] = t09(out)
    sv = sens_values()
    t10(out, sv)
    nums['t10'] = sv
    nums['openkeito_checks'] = extra_numbers()
    vers = sorted({(RUNS / c / 'version.txt').read_text().strip() for c in C.CASES})
    nums['openkeito'] = vers
    (out / 'numbers.json').write_text(json.dumps(nums, indent=1, sort_keys=True) + '\n', encoding='utf-8', newline='\n')
    manifest(out, ', '.join(vers))


def main():
    global RUNS
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--runs', default=str(HERE / 'runs'), help='results of run_cases.py (default: runs)')
    ap.add_argument('--out', default=str(HERE / 'figures'), help='output folder (default: figures)')
    a = ap.parse_args()
    RUNS = Path(a.runs)
    out = Path(a.out)
    build(out)
    print(out, sorted(p.name for p in out.iterdir()))


if __name__ == '__main__':
    main()
