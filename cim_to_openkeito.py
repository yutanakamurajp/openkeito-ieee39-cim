# SPDX-License-Identifier: MIT OR Apache-2.0
# Copyright (c) 2026 The OpenKeito contributors
"""Build OpenKeito models (model format 0.7) of the IEEE 39-bus system directly from CIM.

    python cim_to_openkeito.py <IEEE39.xml> <IEEE39_ic.xml> <out_dir> [--case NAME ...]

Inputs are the CIM instance files of the IEEE 39-bus test system published with emthub / emtiop
(https://github.com/temcdrm/emthub, Apache-2.0): the network (``IEEE39.xml``) and its initial
conditions (``IEEE39_ic.xml``, a power-flow solution with constant-impedance loads).  The CIM is read
with emtiop (``emtiop.api.load_emt_dict``).  Each case of ``CASES`` is written as ``<out_dir>/<case>.json``.

How the network is represented follows the public sources named next to each part:

* Lines (``ACLineSegment``): distributed constant-parameter lines (``bergeron_line``) with the positive-
  and zero-sequence r, x, b of the CIM; a line whose travel time is below the time step would be a
  ``pi_line`` (none in IEEE 39 at 10 us).
* Transformers (``PowerTransformer``): ``three_phase_transformer``.  The CIM tap of end 1 scales the end-1
  rated voltage, and the series impedance (end-1 ohms in the CIM) is given in percent of the nominal end-1
  voltage, so that the impedance seen from end 1 is (r + jx)·tap², as in the CIM/power-flow convention
  ("from voltage / ratio, then impedance").  Magnetizing branches (``TransformerCoreAdmittance``, linear) are
  separate shunt R // L at the end-1 bus, derived as in the star-core function of emtiop (no saturation).
* Loads (``EnergyConsumer``): constant impedance at the nominal voltage (``three_phase_load``), as in the
  emtiop example network export.
* Synchronous machines: IEEE Std 1110 model 2.2 (``synchronous_machine``) with the CIM open-circuit time
  constants.  As in emtiop's machine export (``AppendMachineDynamics``): H = 3.0 s for all
  machines, and Xq'' raised to max(Xl, Xd'') where the CIM value is below Xl or Xd''.  The swing machine
  (bus 31) is a Thevenin source behind Ra + jXd' (emtiop's swing-bus export, ``TheveninVoltage``).
* Excitation and governor of the other machines: the structure of the emtiop example ``syncmach.pch``
  (an ST1A-type static exciter with lead-lag 1+s/1+10s, gain KA from the CIM ``IEEET1``, transient gain
  reduction 1+0.4s/1+0.025s and the limits -4.0·Vt and 4.5·Vt - 0.038·Ifd; a speed governor
  KGOV/(1+0.1s) with KGOV = 1/R and Pmax = max(1.2, VMAX) from the CIM ``TGOV1``; both behind the
  measurement lag 1/(1 + 0.01667 s + 6.944e-5 s²)).
* Photovoltaic plant (bus 30): the WECC generic models REGC_A, REEC_A and REPC_A (WECC Renewable Energy
  Modeling Task Force, "WECC Solar Plant Dynamic Modeling Guidelines" and "Generic Solar Photovoltaic
  System Dynamic Simulation Model Specification", 2014), with the parameters of the CIM dynamics
  (``REGCA1``, ``REECA1``, ``REPCA1``) mapped and initialized as in the emtiop example ``ibr.pch``
  (PLL on the terminal voltage, P/Q references ramped from 0.01 s to 0.02 s, reactive current limit
  enabled from 0.25 s).  The plant is a three-phase controlled current source; an initializing source
  (Thevenin source behind 0.2 pu, as emtiop ``AppendIBRInitializer``) holds the bus until it is opened in
  all three phases at 20 ms.
* PI blocks with limits use ``"anti_windup": "state"`` (the integral itself is held within the limits,
  IEEE Std 421.5 non-windup integrator), the form of the limited PI controllers of ``ibr.pch``.
"""
import argparse
import cmath
import json
import math
from collections import deque
from pathlib import Path

FORMAT_VERSION = '0.7'
FREQ = 60.0
W = 2 * math.pi * FREQ
SYSTEM_ID = '6477751A-0472-4FD6-B3C3-3AD4945CBE56'    # the IEEE 39-bus system in IEEE39.xml
SWING_BUS = '31'
MACHINE_H = 3.0            # emtiop MACHINE_H
IBR_XPU = 0.2              # emtiop AppendIBRInitializer
IBR_T_OPEN = 0.02          # emtiop AppendIBRInitializer (switches open at 20 ms)
IBR_THETA0 = 1.8           # emtiop AppendSolar (initial PLL angle)
MIN_IMAG, MIN_PNLL = 0.00015, 0.00010   # emtiop (core without b or g)
# syncmach.pch exciter and governor constants (emtiop examples)
TC, TB, TLED, TLAG = 1.0, 10.0, 0.4, 0.025
VRMIN, VRMAX, KC = 4.0, 4.5, 0.038
EXCLIM_MAX = 1.5
T_GOV = 0.1                # KGOV (1 + 5 s)/(1 + 5.1 s + 0.5 s^2) = KGOV/(1 + 0.1 s)
LAG2 = (0.01667, 6.944E-5)  # measurement lag 1/(1 + a1 s + a2 s^2) of syncmach.pch and ibr.pch
PLL_KP, PLL_KI = 188.5, 1.34E5   # ibr.pch PLL
WC = 376.9911              # ibr.pch nominal angular frequency
# ibr.pch lookup tables (REGC_A LVPL and low-voltage active current management, REEC_A voltage-dip
# reactive current, VDL1/VDL2, REPC_A error limiter)
TABLES = {'LVPL': [(-1, 0), (0.4, 0), (0.9, 1.2), (10, 1.2)],
          'LVACM': [(-1, 0), (0.4, 0), (0.8, 1), (10, 1)],
          'IQVRF': [(-1e8, -1.05), (-1.05, -1.05), (-0.05, -0.05), (-0.0499, 0), (0.0499, 0), (0.05, 0.05),
                    (1.05, 1.05), (1e8, 1.05)],
          'VDL': [(-1, 1), (2, 1), (3, 1), (4, 1)],
          'QEX2': [(-2000, -999), (-999, -999), (999, 999), (2000, 999)]}


# ---------------------------------------------------------------- CIM

def load_cim(cim_path, ic_path, system_id=SYSTEM_ID):
    """Read the CIM network and its initial conditions with emtiop into one table dictionary."""
    import rdflib
    import emtiop.api as emtiop
    g = rdflib.Graph(store='Oxigraph')
    g.parse(str(cim_path))
    g.parse(str(ic_path))
    return emtiop.load_emt_dict(g, system_id)


def network(d):
    """Collect buses, lines, transformers, machines, PV plants and loads (SI units) and the initial
    conditions (bus voltage magnitude in pu and angle in degrees; machine and PV power in MW and Mvar)."""
    buses = {k: {'name': v['name'], 'kv': v['BaseVoltage_nominalVoltage'] / 1e3} for k, v in d['EMTBus']['vals'].items()}
    name = lambda k: buses[k]['name']
    lines = [{'name': v['name'], 'from': name(v['ConnectivityNode1_mRID']), 'to': name(v['ConnectivityNode2_mRID']),
              'length_m': v['length'], 'r': v['r'], 'x': v['x'], 'bch': v['bch'], 'r0': v['r0'], 'x0': v['x0'],
              'b0ch': v['b0ch']} for v in d['EMTLine']['vals'].values()]
    ends = {}
    for k, v in d['EMTPowerXfmrWinding']['vals'].items():
        pid, n = k.split(':')
        ends.setdefault(pid, {})[int(n)] = v
    xfmrs = []
    for pid, e in ends.items():
        if len(e) != 2:
            raise NotImplementedError(f'{pid}: only two-winding transformers')
        e1, e2 = e[1], e[2]
        mesh = d['EMTPowerXfmrMesh']['vals'][f"{e1['mRID']}:{e2['mRID']}"]
        tap = 1.0
        t = d['EMTXfmrTap']['vals'].get(e1['mRID'])
        if t:
            tap = 1.0 + t['stepVoltageIncrement'] * t['step'] / 100.0
        if d['EMTXfmrTap']['vals'].get(e2['mRID']):
            raise NotImplementedError(f'{pid}: tap on end 2')
        core = d.get('EMTPowerXfmrCore', {}).get('vals', {}).get(pid)
        xfmrs.append({'name': e1['PowerTransformer_name'], 'from': name(e1['ConnectivityNode1_mRID']),
                      'to': name(e2['ConnectivityNode1_mRID']), 'kv1': e1['ratedU'] / 1e3, 'kv2': e2['ratedU'] / 1e3,
                      'mva': e1['ratedS'] / 1e6, 'r': mesh['r'], 'x': mesh['x'], 'tap': tap,
                      'conn1': e1['connectionKind'], 'conn2': e2['connectionKind'],
                      'clock': e2['phaseAngleClock'] - e1['phaseAngleClock'],
                      'grounded1': str(e1.get('grounded')).lower() == 'true',
                      'grounded2': str(e2.get('grounded')).lower() == 'true',
                      'core': core and {'b': core['b'], 'g': core['g'],
                                        'kv': e[core['endNumber']]['ratedU'] / 1e3}})
    types = {k: v['name'].upper() for k, v in d.get('EMTDynamicsModelType', {}).get('vals', {}).items()}
    dyn = {}
    for key, row in d.get('EMTDynamicsParameter', {}).get('vals', {}).items():
        eq, tid, prm = key.split(':')
        dyn.setdefault(eq, {}).setdefault(types.get(tid, tid), {})[prm] = row['value']
    bv, bf = d['EMTBusVoltage']['vals'], d['EMTBranchFlow']['vals']
    ic_bus = {name(k): (bv[k]['v'] / (b['kv'] * 1e3), bv[k]['angle']) for k, b in buses.items()}
    gens = []
    for k, v in d['EMTSyncMachine']['vals'].items():
        f = bf[f'{k}:1']                     # terminal flow in load convention
        gens.append(dict(v, bus=name(v['ConnectivityNode1_mRID']), dyn=dyn.get(k, {}),
                         p_mw=-f['p'] / 1e6, q_mvar=-f['q'] / 1e6))
    pvs = []
    for k, v in d.get('EMTSolar', {}).get('vals', {}).items():
        f = bf[f'{k}:1']
        pvs.append(dict(v, bus=name(v['ConnectivityNode1_mRID']), dyn=dyn.get(k, {}),
                        p_mw=-f['p'] / 1e6, q_mvar=-f['q'] / 1e6))
    if d.get('EMTWind', {}).get('vals'):
        raise NotImplementedError('wind plants')
    loads = [{'name': v['name'], 'bus': name(v['ConnectivityNode1_mRID']), 'kv': v['BaseVoltage_nominalVoltage'] / 1e3,
              'p_mw': v['p'] / 1e6, 'q_mvar': v['q'] / 1e6} for v in d['EMTLoad']['vals'].values()]
    net = {'buses': {b['name']: b for b in buses.values()}, 'lines': lines, 'xfmrs': xfmrs, 'gens': gens,
           'pvs': pvs, 'loads': loads, 'ic': ic_bus}
    net['shift'] = phase_shift(net, SWING_BUS)
    return net


def phase_shift(net, swing):
    """Angle shift of every bus against the swing bus from the transformer clock numbers (the end-2 side
    lags by clock x 30 deg).  The CIM power flow does not contain these shifts; the three-phase model does."""
    shift = {swing: 0.0}
    q = deque([swing])
    while q:
        u = q.popleft()
        for br in net['xfmrs'] + net['lines']:
            s = -30.0 * br.get('clock', 0)
            for a, b, sign in ((br['from'], br['to'], 1), (br['to'], br['from'], -1)):
                if a == u and b not in shift:
                    shift[b] = shift[u] + sign * s
                    q.append(b)
    return shift


# ---------------------------------------------------------------- network components

def bus_id(name):
    return f'B{name}'


def vector_group(t):
    hv = {'Y': 'YN' if t['grounded1'] else 'Y', 'D': 'D'}[t['conn1']]
    lv = {'Y': 'yn' if t['grounded2'] else 'y', 'D': 'd'}[t['conn2']]
    return f"{hv}{lv}{t['clock'] % 12}"


def line(ln, dt):
    km = ln['length_m'] / 1e3
    per_km = {'type': 'sequence',
              'r1_ohm_per_km': ln['r'] / km, 'l1_h_per_km': ln['x'] / W / km, 'c1_f_per_km': ln['bch'] / W / km,
              'r0_ohm_per_km': ln['r0'] / km, 'l0_h_per_km': ln['x0'] / W / km, 'c0_f_per_km': ln['b0ch'] / W / km}
    tau = min(km * math.sqrt(per_km['l1_h_per_km'] * per_km['c1_f_per_km']),
              km * math.sqrt(per_km['l0_h_per_km'] * per_km['c0_f_per_km']))
    c = {'type': 'bergeron_line', 'name': f"L{ln['name']}", 'buses': [bus_id(ln['from']), bus_id(ln['to'])],
         'length_km': km, 'per_km': per_km}
    if tau < dt:
        c.update(type='pi_line', sections=1)
    return c, tau


def transformer(t):
    zb1 = t['kv1'] ** 2 / t['mva']
    buses, kv = [bus_id(t['from']), bus_id(t['to'])], [t['kv1'] * t['tap'], t['kv2']]
    if t['kv1'] < t['kv2']:
        # OpenKeito takes the high-voltage side first; a Y-Y (clock 0) bank is symmetric, so swap the sides
        if not (t['conn1'] == t['conn2'] == 'Y' and t['clock'] % 12 == 0):
            raise NotImplementedError(f"{t['name']}: end 1 on the low-voltage side")
        buses, kv = buses[::-1], kv[::-1]
    return {'type': 'three_phase_transformer', 'name': f"T{t['name']}", 'buses': buses, 'vector_group': vector_group(t),
            'rated_voltage_kv': kv, 'rated_power_mva': t['mva'], 'rated_frequency_hz': FREQ,
            'impedance_percent': 100 * math.hypot(t['r'], t['x']) / zb1, 'resistance_percent': 100 * t['r'] / zb1}


def magnetizing(t):
    """Linear magnetizing branch of end 1 as R and L per phase to ground (emtiop's star-core function without
    saturation: the CIM b, g of the core end referred to end 1; Imag = sqrt2·b·Vph, Fmag = sqrt2·Vph/omega)."""
    core = t['core']
    if t['conn1'] != 'Y':
        raise NotImplementedError(f"{t['name']}: magnetizing branch of a delta end 1")
    ratio = (t['kv1'] / core['kv']) ** 2
    b, g = core['b'] / ratio, core['g'] / ratio
    v_ph, s_ph = t['kv1'] * 1e3 / math.sqrt(3.0), t['mva'] * 1e6 / 3.0
    imag = math.sqrt(2.0) * b * v_ph if b > 0 else MIN_IMAG * math.sqrt(2.0) * s_ph / v_ph
    fmag = v_ph / W * math.sqrt(2.0)
    rmag = 1.0 / g if g > 0 else v_ph ** 2 / MIN_PNLL / s_ph
    return rmag, fmag / imag


def magnetizing_percent(t):
    """The same branch as no-load current and loss of a three_phase_transformer (rated voltage kv1·tap)."""
    r, l = magnetizing(t)
    v0 = t['kv1'] * t['tap'] * 1e3
    return 100 * v0 ** 2 * math.hypot(1 / r, 1 / (W * l)) / (t['mva'] * 1e6), v0 ** 2 / r / 1e3


def load(ld):
    return {'type': 'three_phase_load', 'name': f"LD{ld['name']}", 'bus': bus_id(ld['bus']), 'rated_voltage_kv': ld['kv'],
            'rated_frequency_hz': FREQ, 'active_power_mw': ld['p_mw'], 'reactive_power_mvar': ld['q_mvar']}


def xq2(g):
    """Xq'' as in emtiop AppendMachineDynamics."""
    xl, xd2, xq2_ = g['statorLeakageReactance'], g['xDirectSubtrans'], g['xQuadSubtrans']
    if xq2_ < xl or xq2_ < xd2:
        xq2_ = max(xl, xd2)
    if not xl < min(xd2, xq2_):
        raise ValueError(f"{g['name']}: Xl must be below Xd'' and Xq''")
    return xq2_


def machine(g, net):
    v, a = net['ic'][g['bus']]
    kv_bus = net['buses'][g['bus']]['kv']
    rated_kv = g['ratedU'] / 1e3
    return {'type': 'synchronous_machine', 'name': f"G{g['name']}", 'bus': bus_id(g['bus']),
            'rated_power_mva': g['ratedS'] / 1e6, 'rated_voltage_kv': rated_kv, 'frequency_hz': FREQ,
            'xd_pu': g['xDirectSync'], 'xq_pu': g['xQuadSync'], 'xd_transient_pu': g['xDirectTrans'],
            'xq_transient_pu': g['xQuadTrans'], 'xd_subtransient_pu': g['xDirectSubtrans'], 'xq_subtransient_pu': xq2(g),
            'xl_pu': g['statorLeakageReactance'], 'ra_pu': g['statorResistance'],
            'td0_transient_s': g['tpdo'], 'td0_subtransient_s': g['tppdo'],
            'tq0_transient_s': g['tpqo'], 'tq0_subtransient_s': g['tppqo'],
            'inertia_h_s': MACHINE_H, 'damping_pu': 0.0,
            'operating_point': {'active_power_mw': g['p_mw'], 'reactive_power_mvar': g['q_mvar'],
                                'voltage_pu': v * kv_bus / rated_kv, 'angle_deg': a + net['shift'][g['bus']]}}


def thevenin(vpu, deg, rpu, xpu, ppu, qpu):
    """Source voltage behind r + jx that delivers p + jq at the terminal voltage (emtiop TheveninVoltage)."""
    vt = cmath.rect(vpu, math.radians(deg))
    vs = vt + complex(rpu, xpu) * (complex(ppu, qpu) / vt).conjugate()
    return abs(vs), math.degrees(cmath.phase(vs))


def swing_source(g, net):
    v, a = net['ic'][g['bus']]
    s = g['ratedS'] / 1e6
    kv = net['buses'][g['bus']]['kv']
    zb = kv ** 2 / s
    e, ang = thevenin(v, a, g['statorResistance'], g['xDirectTrans'], g['p_mw'] / s, g['q_mvar'] / s)
    r, l = g['statorResistance'] * zb, g['xDirectTrans'] * zb / W
    return {'type': 'three_phase_source', 'name': f"SW{g['name']}", 'bus': bus_id(g['bus']),
            'line_voltage_rms_v': e * kv * 1e3, 'frequency_hz': FREQ, 'angle_deg': ang + net['shift'][g['bus']],
            'impedance': {'type': 'sequence', 'r1_ohm': r, 'l1_h': l, 'r0_ohm': r, 'l0_h': l}}


def pv_initializer(p, net, bus):
    v, a = net['ic'][p['bus']]
    s = p['ratedS'] / 1e6
    kv = net['buses'][p['bus']]['kv']
    e, ang = thevenin(v, a, 0.0, IBR_XPU, p['p_mw'] / s, p['q_mvar'] / s)
    l = IBR_XPU * kv ** 2 / s / W
    return {'type': 'three_phase_source', 'name': f"IBRS{p['name']}", 'bus': bus,
            'line_voltage_rms_v': e * kv * 1e3, 'frequency_hz': FREQ, 'angle_deg': ang + net['shift'][p['bus']],
            'impedance': {'type': 'sequence', 'r1_ohm': 0.0, 'l1_h': l, 'r0_ohm': 0.0, 'l0_h': l}}


# ---------------------------------------------------------------- control blocks

class Blocks:
    """Control blocks whose names start with a common prefix."""

    def __init__(self, prefix):
        self.p, self.items = prefix, []

    def n(self, name):
        return f'{self.p}{name}'

    def ref(self, x):
        return self.n(x) if isinstance(x, str) else float(x)

    def add(self, typ, name, **kw):
        self.items.append(dict({'type': typ, 'name': self.n(name)}, **kw))
        return self.n(name)

    def const(self, name, v):
        return self.add('constant', name, value=float(v))

    def gain(self, name, src, k):
        return self.add('gain', name, input=self.ref(src), gain=float(k))

    def sum(self, name, *terms):
        return self.add('sum', name, inputs=[self.ref(s) for s, _ in terms], signs=[float(k) for _, k in terms])

    def lin(self, name, *terms):
        """sum of k·signal (a gain block where k is not +-1)."""
        ins = []
        for i, (s, k) in enumerate(terms):
            ins.append((s, k) if abs(k) == 1.0 else (self.gain(f'{name}_k{i}', s, k)[len(self.p):], 1.0))
        return self.sum(name, *ins)

    def prod(self, name, *src):
        return self.add('product', name, inputs=[self.ref(s) for s in src])

    def quot(self, name, num, den, minimum=None):
        kw = {'minimum_denominator': minimum} if minimum else {}
        return self.add('quotient', name, inputs=[self.ref(num), self.ref(den)], **kw)

    def fn(self, name, f, *src):
        return self.add('function', name, function=f, inputs=[self.ref(s) for s in src])

    def lag(self, name, src, t, init=None, lo=None, hi=None, k=1.0):
        kw = {}
        if init is not None:
            kw['initial_input'] = float(init)
        if lo is not None:
            kw['minimum'] = self.ref(lo)
        if hi is not None:
            kw['maximum'] = self.ref(hi)
        if k != 1.0:
            kw['gain'] = float(k)
        return self.add('lag', name, input=self.ref(src), time_constant_s=float(t), **kw)

    def pi(self, name, err, kp, ki, init=0.0, lo=None, hi=None, anti_windup=None):
        kw = {'anti_windup': anti_windup} if anti_windup else {}
        if lo is not None:
            kw['output_min'] = self.ref(lo)
        if hi is not None:
            kw['output_max'] = self.ref(hi)
        return self.add('pi', name, input=self.ref(err), proportional_gain=float(kp), integral_gain_per_s=float(ki),
                        initial_integral=float(init), **kw)

    def limiter(self, name, src, lo=None, hi=None):
        kw = {}
        if lo is not None:
            kw['minimum'] = self.ref(lo)
        if hi is not None:
            kw['maximum'] = self.ref(hi)
        return self.add('limiter', name, input=self.ref(src), **kw)

    def table(self, name, src, pts):
        return self.add('piecewise_linear', name, input=self.ref(src), x=[float(p[0]) for p in pts],
                        y=[float(p[1]) for p in pts])

    def meas(self, name, comp, signal):
        return self.add('measurement', name, quantity={'type': 'component_signal', 'component': comp, 'signal': signal})


def lag2_taus():
    a1, a2 = LAG2
    r = math.sqrt(a1 ** 2 - 4 * a2)
    return (a1 + r) / 2, (a1 - r) / 2


def machine_parameters(g):
    """KGOV, PMAX and KA as in emtiop AppendMachineDynamics (TGOV1 and IEEET1 in the CIM)."""
    dyn = g.get('dyn', {})
    kgov, pmax = 25.0, 1.2
    if 'TGOV1' in dyn:
        kgov, pmax = 1.0 / dyn['TGOV1']['R'], max(1.2, dyn['TGOV1']['VMAX'])
    if 'IEEET1' not in dyn:
        raise ValueError(f"{g['name']}: no IEEET1 (exciter gain)")
    return kgov, pmax, dyn['IEEET1']['KA']


def machine_controls(g, dt):
    """Exciter and governor of one machine (syncmach.pch structure).  Returns the blocks and the signals of
    the field voltage and the mechanical torque.  The references hold the initial operating point."""
    gid = f"G{g['name']}"
    b = Blocks(f'{gid}_')
    t1, t2 = lag2_taus()
    kgov, pmax, ka = machine_parameters(g)
    b.meas('EA', gid, 'terminal_voltage_pu')
    b.meas('OMG', gid, 'speed_pu')
    b.meas('IFD', gid, 'field_current_pu')
    b.meas('EFS', gid, 'initial_field_voltage_pu')
    b.meas('TGI', gid, 'initial_mechanical_torque_pu')
    b.add('initial_value', 'EAS', input=b.n('EA'))
    b.add('initial_value', 'IFD0', input=b.n('IFD'))
    b.const('ONE', 1.0)
    # exciter
    b.lag('VRMS1', 'EA', t1)
    b.lag('VRMS', 'VRMS1', t2)
    b.const('INVKA', 1.0 / ka)
    b.sum('EXL0', ('EAS', 1), ('INVKA', 1), ('VRMS', -1))           # V0 = Vt(0) + 1/KA
    b.limiter('EXCLIM', 'EXL0', -EXCLIM_MAX, EXCLIM_MAX)
    b.lag('EXTL', 'EXCLIM', TB)                                     # (1 + TC s)/(1 + TB s)
    b.lin('EXCVA', ('EXCLIM', ka * TC / TB), ('EXTL', ka * (1 - TC / TB)))
    b.lag('EXOL', 'EXCVA', TLAG)                                    # (1 + TLED s)/(1 + TLAG s)
    b.lin('EXCOUT', ('EXCVA', TLED / TLAG), ('EXOL', 1 - TLED / TLAG))
    b.quot('IFDN0', 'IFD', 'IFD0')                                  # field current in units of its initial value
    b.add('delay', 'IFDN', input=b.n('IFDN0'), delay_s=float(dt), initial_output=1.0)
    b.lin('EFU', ('VRMS', VRMAX), ('IFDN', -KC))
    b.gain('EFL', 'VRMS', -VRMIN)
    b.fn('EFD1', 'max', 'EFL', 'EXCOUT')
    b.fn('EFD', 'min', 'EFU', 'EFD1')
    b.prod('EFO', 'EFD', 'EFS')
    # governor
    b.sum('DW0', ('ONE', 1), ('OMG', -1))
    b.lag('DW1', 'DW0', t1)
    b.lag('DW', 'DW1', t2)
    b.lag('DPWG', 'DW', T_GOV, k=kgov)
    b.sum('PM0', ('ONE', 1), ('DPWG', 1))
    b.const('PMAX', pmax)
    b.fn('PMECH', 'min', 'PMAX', 'PM0')
    b.prod('TGP', 'TGI', 'PMECH')
    b.quot('TMO', 'TGP', 'OMG')
    return b.items, b.n('EFO'), b.n('TMO')


def pv_parameters(p):
    """Parameters of REGC_A / REEC_A / REPC_A from the CIM dynamics, mapped as in emtiop AppendSolar."""
    reec, regc, repc = p['dyn']['REECA1'], p['dyn']['REGCA1'], p['dyn']['REPCA1']
    s, v = p['ratedS'], p['ratedU']
    ppu, qpu = p['p_mw'] * 1e6 / s, p['q_mvar'] * 1e6 / s
    return {'VOLTBASE': v, 'POWERBASE': s, 'CURRBASE': s / v, 'PPU': ppu, 'QPU': qpu, 'VREFPU': 1.0,
            'ICTHETA': IBR_THETA0, 'VFLAG': int(reec['VFLAG']), 'QFLAG': int(reec['QFLAG']), 'PQFLG': int(reec['PQFLAG']),
            'PFFLG': int(reec['PFFLAG']), 'REFFG': int(repc['RefFlag']), 'LDCFG': int(repc['VCFlag']),
            'FRQFG': int(repc['Fflag']), 'LVPSW': int(regc['Lvplsw']), 'TANPF': qpu / ppu if ppu > 0 else 0.0,
            'PFLAG': int(reec['PFLAG']), 'PMIN': reec['PMIN'], 'PMAX': reec['PMAX'], 'IQRMX': regc['Iqrmax'],
            'IQRMN': regc['Iqrmin'], 'VDIP': reec['Vdip'], 'VUP': reec['Vup'], 'KQV': reec['Kqv'],
            'KC': repc['Kc'], 'ZC2': repc['Rc'] ** 2 + repc['Xc'] ** 2, 'VQMAX': reec['VMAX'], 'VQMIN': reec['VMIN'],
            'IMAX': reec['Imax'], 'QMIN': repc['Qmin'], 'QMAX': repc['Qmax'], 'KQP': max(1e-8, reec['Kqp']),
            'KQI': reec['Kqi'], 'KVP': max(1e-8, reec['Kvp']), 'KVI': reec['Kvi'], 'KP': max(1e-8, repc['Kp']),
            'KI': repc['Ki']}


def pv_controls(a, src, bus, qreg0, pi_form='state'):
    """REGC_A, REEC_A and REPC_A of one PV plant as control blocks (the equations of ibr.pch).  Terms that
    the flags of this plant make zero are left out (other flag values are refused).  pi_form: 'state'
    (default, limited integral as ibr.pch) or 'conditional' (OpenKeito's default PI, for a sensitivity)."""
    want = {'PFFLG': 0, 'PFLAG': 0, 'PQFLG': 0, 'QFLAG': 1, 'VFLAG': 1, 'REFFG': 1, 'LDCFG': 1, 'LVPSW': 1, 'KC': 0.0,
            'ZC2': 0.0}
    bad = {k: a[k] for k, v in want.items() if a[k] != v}
    if bad:
        raise NotImplementedError(f'flags not covered by this conversion: {bad}')
    b = Blocks(f'{src}_')
    aw = 'state' if pi_form == 'state' else None
    t1, t2 = lag2_taus()
    vb, ib = a['VOLTBASE'], a['CURRBASE']
    b.add('three_phase_measurement', 'V', quantity={'type': 'bus_voltage', 'bus': bus})
    b.add('three_phase_measurement', 'I', quantity={'type': 'branch_current', 'component': src})
    # PLL: THETA = integral of WE (wrapped to +-pi), WE = WC - PI(VGQ), starting at ICTHETA
    b.add('integrator', 'THETA', input=b.n('WE'), initial_value=float(a['ICTHETA']), minimum=-math.pi, maximum=math.pi,
          wrap=True)
    b.add('abc_to_dq0', 'VDQ', inputs=[b.n('V.a'), b.n('V.b'), b.n('V.c')], angle=b.n('THETA'))
    b.add('abc_to_dq0', 'IDQ', inputs=[b.n('I.a'), b.n('I.b'), b.n('I.c')], angle=b.n('THETA'))
    # ibr.pch's alpha-beta scaling: VGD = sqrt(3/2)·d / VOLTBASE (d of the amplitude-invariant transform)
    b.gain('VGD', 'VDQ.d', math.sqrt(1.5) / vb)
    b.gain('VGQ', 'VDQ.q', math.sqrt(1.5) / vb)
    b.gain('IGD', 'IDQ.d', math.sqrt(1.5) / ib)
    b.gain('IGQ', 'IDQ.q', math.sqrt(1.5) / ib)
    b.pi('DW', 'VGQ', PLL_KP, PLL_KI)
    b.const('WC', WC)
    b.sum('WE', ('WC', 1), ('DW', -1))
    b.prod('PD1', 'VGD', 'IGD')
    b.prod('PD2', 'VGQ', 'IGQ')
    b.sum('PDQ0', ('PD1', 1), ('PD2', 1))
    b.prod('QD1', 'VGQ', 'IGD')
    b.prod('QD2', 'VGD', 'IGQ')
    b.sum('QDQ0', ('QD1', 1), ('QD2', -1))
    b.fn('VGR', 'hypot', 'VGD', 'VGQ')
    b.lag('VGPU1', 'VGR', t1, init=a['VREFPU'])
    b.lag('VGPU', 'VGPU1', t2, init=a['VREFPU'])
    b.lag('PG1', 'PDQ0', t1, init=0.0)
    b.lag('PG', 'PG1', t2, init=0.0)
    b.lag('QG1', 'QDQ0', t1, init=0.0)
    b.lag('QG', 'QG1', t2, init=0.0)
    # REEC_A inputs
    b.lag('VFILT', 'VGPU', 0.03, init=a['VREFPU'])
    b.const('VREF0', a['VREFPU'])
    b.const('VDIP', a['VDIP'])
    b.const('VUP', a['VUP'])
    b.add('comparator', 'VLO', inputs=[b.n('VDIP'), b.n('VGPU')])
    b.add('comparator', 'VHI', inputs=[b.n('VGPU'), b.n('VUP')])
    b.sum('VOUT', ('VLO', 1), ('VHI', 1))
    b.const('ONE', 1.0)
    b.fn('VOUT1', 'min', 'ONE', 'VOUT')
    b.gain('VSUP', 'VOUT1', a['KQV'])
    b.sum('IQVRFin', ('VREF0', 1), ('VFILT', -1))
    b.table('IQVRF', 'IQVRFin', TABLES['IQVRF'])
    b.table('VDL1', 'VFILT', TABLES['VDL'])
    b.table('VDL2', 'VFILT', TABLES['VDL'])
    b.add('time_series', 'RAMP', times_s=[0.0, 0.01, 0.02], values=[0.0, 0.0, 1.0], interpolation='linear')
    b.gain('PREF', 'RAMP', a['PPU'])
    # REPC_A (voltage control with line-drop compensation off: VCTRL = VFILT)
    b.lag('VCLAG', 'VFILT', 0.02, init=a['VREFPU'])
    b.sum('QEX1', ('VREF0', 1), ('VCLAG', -1))
    b.table('QEX2', 'QEX1', TABLES['QEX2'])
    b.pi('QEX3', 'QEX2', a['KP'], a['KI'], init=a['QPU'], lo=a['QMIN'], hi=a['QMAX'], anti_windup=aw)
    b.lag('QEX4', 'QEX3', 0.08, init=a['QPU'])
    # REEC_A reactive and current commands
    b.limiter('QCTRL', 'QEX4', a['QMIN'], a['QMAX'])
    b.sum('QERR', ('QCTRL', 1), ('QG', -1))
    b.pi('QREG', 'QERR', a['KQP'], a['KQI'], init=qreg0, lo=a['VQMIN'], hi=a['VQMAX'], anti_windup=aw)
    b.add('step', 'IQLIM', time_s=0.25, before=0.0, after=float(a['IMAX']))
    b.fn('IQMAX', 'min', 'VDL1', 'IQLIM')
    b.gain('IQMIN', 'IQMAX', -1.0)
    b.sum('VERR', ('QREG', 1), ('VFILT', -1))
    b.pi('IQRAW', 'VERR', a['KVP'], a['KVI'], init=0.0, lo='IQMIN', hi='IQMAX', anti_windup=aw)
    b.prod('IQINJ', 'IQVRF', 'VSUP')
    b.sum('IQSUM', ('IQINJ', 1), ('IQRAW', 1))
    b.limiter('IQCMD', 'IQSUM', 'IQMIN', 'IQMAX')
    b.prod('IQ2', 'IQCMD', 'IQCMD')
    b.const('IMAX2', a['IMAX'] ** 2)
    b.sum('IPL2', ('IMAX2', 1), ('IQ2', -1))
    b.const('ZERO', 0.0)
    b.fn('IPL2p', 'max', 'ZERO', 'IPL2')
    b.fn('IPLIM', 'sqrt', 'IPL2p')
    b.fn('IPMAX', 'min', 'VDL2', 'IPLIM')
    b.lag('PORD', 'PREF', 0.03, init=0.0, lo=a['PMIN'], hi=a['PMAX'])
    b.quot('IPRAW', 'PORD', 'VFILT', minimum=1e-3)
    b.limiter('IPCMD', 'IPRAW', a['PMIN'], 'IPMAX')
    # REGC_A
    b.lag('VLVPL', 'VGPU', 0.015, init=0.0)
    b.table('LVPL', 'VLVPL', TABLES['LVPL'])
    b.table('LVACM', 'VGPU', TABLES['LVACM'])
    b.lag('IPREF', 'IPCMD', 0.02, init=0.0, lo=-1e9, hi='LVPL')
    b.lag('IQREF', 'IQCMD', 0.02, init=0.0, lo=a['IQRMN'], hi=a['IQRMX'])
    b.prod('IDL', 'IPREF', 'LVACM')
    b.gain('ID', 'IDL', -1.0)
    b.add('dq0_to_abc', 'JDQ', d=b.n('ID'), q=b.n('IQREF'), angle=b.n('THETA'))
    sig = [b.gain(f'J{ph.upper()}', f'JDQ.{ph}', math.sqrt(2 / 3) * ib) for ph in 'abc']
    return b.items, sig, b


RECORD_PV = ('PDQ0', 'QDQ0', 'WE', 'VGPU', 'VFILT', 'QREG', 'QEX3', 'QEX4', 'IQCMD', 'IPCMD', 'IQRAW', 'PG', 'QG',
             'VSUP', 'LVPL', 'IPREF', 'IQREF')
RECORD_GEN = ('OMG', 'EA', 'TMO', 'EFO')


# ---------------------------------------------------------------- cases

CASES = {
    # no fault, 6 s; the REEC_A voltage PI (QREG) starts at the bus voltage, or at 0 as in ibr.pch (then held at VMIN)
    'nofault': {'end_time_s': 6.0, 'output_interval_steps': 50, 'qreg0': 'bus_voltage', 'magnetizing': None, 'fault': None,
              'record_currents': True},
    'nofault_qreg0': {'end_time_s': 6.0, 'output_interval_steps': 50, 'qreg0': 'zero', 'magnetizing': None, 'fault': None,
               'record_currents': True},
    # single line to ground, bus 16 phase C, 2.0-2.2 s
    'slg16': {'end_time_s': 3.5, 'output_interval_steps': 10, 'qreg0': 'bus_voltage', 'magnetizing': None,
           'fault': {'bus': '16', 'phases': 'c', 'r_ohm': 0.001, 't_on': 2.0, 't_off': 2.2}},
    # three-phase fault, bus 2, 3.0-3.1 s
    'tpf02': {'end_time_s': 4.0, 'output_interval_steps': 50, 'qreg0': 'bus_voltage', 'magnetizing': None,
                'fault': {'bus': '2', 'phases': 'abc', 'r_ohm': 0.001, 't_on': 3.0, 't_off': 3.1}},
    # 20 s, the fault of the emtiop example IEEE39.prm (bus 16 phase C, 12.0-12.2 s), magnetizing branches
    # connected from t = 0+ (excluded from the phasor solution)
    'slg16_20s': {'end_time_s': 20.0, 'output_interval_steps': 50, 'qreg0': 'zero', 'magnetizing': 'shunt_exclude',
                   'fault': {'bus': '16', 'phases': 'c', 'r_ohm': 0.001, 't_on': 12.0, 't_off': 12.2}},
    # sensitivities: OpenKeito's conditional-integration PI; magnetizing inside the transformers or as
    # shunts in the phasor solution
    'nofault_conditional_pi': {'end_time_s': 6.0, 'output_interval_steps': 50, 'qreg0': 'bus_voltage', 'magnetizing': None, 'fault': None,
                   'pi_form': 'conditional'},
    'slg16_20s_mag_transformer': {'end_time_s': 20.0, 'output_interval_steps': 50, 'qreg0': 'zero', 'magnetizing': 'transformer',
                 'fault': {'bus': '16', 'phases': 'c', 'r_ohm': 0.001, 't_on': 12.0, 't_off': 12.2}},
    'slg16_20s_mag_shunt': {'end_time_s': 20.0, 'output_interval_steps': 50, 'qreg0': 'zero', 'magnetizing': 'shunt',
                  'fault': {'bus': '16', 'phases': 'c', 'r_ohm': 0.001, 't_on': 12.0, 't_off': 12.2}},
}


def build_model(net, case, name, dt=1e-5, method='2s-dirk'):
    """The OpenKeito model (dict) of one case."""
    comps, controls, outputs = [], [], []
    for ln in net['lines']:
        comps.append(line(ln, dt)[0])
    mag = case.get('magnetizing')
    for t in net['xfmrs']:
        c = transformer(t)
        comps.append(c)
        if mag == 'transformer':
            c['no_load_current_percent'], c['no_load_loss_kw'] = magnetizing_percent(t)
        elif mag in ('shunt', 'shunt_exclude'):
            r, l = magnetizing(t)
            kv = t['kv1']
            m = {'type': 'three_phase_load', 'name': f"MAG{t['name']}", 'bus': bus_id(t['from']), 'rated_voltage_kv': kv,
                 'rated_frequency_hz': FREQ, 'active_power_mw': kv ** 2 / r, 'reactive_power_mvar': kv ** 2 / (W * l)}
            if mag == 'shunt_exclude':
                m['initialization'] = 'exclude'
            comps.append(m)
        elif mag is not None:
            raise ValueError(f'magnetizing: {mag}')
    comps += [load(ld) for ld in net['loads']]
    currents = []
    for g in net['gens']:
        if g['bus'] == SWING_BUS:
            sw = swing_source(g, net)
            comps.append(sw)
            currents.append(sw['name'])
            b = Blocks('SWP_')        # power of the swing source (bus to source current: minus the delivered power)
            b.add('three_phase_measurement', 'V', quantity={'type': 'bus_voltage', 'bus': sw['bus']})
            b.add('three_phase_measurement', 'I', quantity={'type': 'branch_current', 'component': sw['name']})
            for ph in 'abc':
                b.prod(f'P{ph}', f'V.{ph}', f'I.{ph}')
            b.lin('W', ('Pa', 1e-6), ('Pb', 1e-6), ('Pc', 1e-6))
            controls += b.items
            outputs.append({'type': 'signal', 'signal': 'SWP_W'})
            continue
        m = machine(g, net)
        items, efd, tm = machine_controls(g, dt)
        m['field_voltage_signal'], m['mechanical_torque_signal'] = efd, tm
        comps.append(m)
        controls += items
        currents.append(m['name'])
        outputs += [{'type': 'signal', 'signal': f"{m['name']}_{s}"} for s in RECORD_GEN]
        controls.append({'type': 'measurement', 'name': f"{m['name']}_P",
                         'quantity': {'type': 'component_signal', 'component': m['name'], 'signal': 'active_power_pu'}})
        outputs.append({'type': 'signal', 'signal': f"{m['name']}_P"})
    for p in net['pvs']:
        bus = bus_id(p['bus'])
        sbus = f"S{p['bus']}"
        init = pv_initializer(p, net, sbus)
        comps.append(init)
        currents.append(init['name'])
        for ph in 'abc':
            comps.append({'type': 'switch', 'name': f"{init['name']}_SW{ph}", 'nodes': [f'{sbus}.{ph}', f'{bus}.{ph}'],
                          'initially_closed': True, 'open_time_s': IBR_T_OPEN, 'open_at_current_zero': False})
        cid = f"IBR{p['name']}"
        a = pv_parameters(p)
        qreg0 = net['ic'][p['bus']][0] if case['qreg0'] == 'bus_voltage' else 0.0
        items, sig, blk = pv_controls(a, cid, bus, qreg0, case.get('pi_form', 'state'))
        comps.append({'type': 'three_phase_controlled_current_source', 'name': cid, 'bus': bus, 'signals': sig,
                      'steady_state': {'current_rms_a': 0.0, 'angle_deg': 0.0}})
        controls += items
        outputs += [{'type': 'signal', 'signal': blk.n(s)} for s in RECORD_PV]
    f = case.get('fault')
    if f:
        fid = f"F{f['bus']}"
        comps.append({'type': 'fault', 'name': fid, 'bus': bus_id(f['bus']), 'phases': list(f['phases']), 'ground': True,
                      'fault_resistance_ohm': f['r_ohm'], 'ground_resistance_ohm': 0.0, 'start_time_s': f['t_on'],
                      'clear_time_s': f['t_off'], 'clear_at_current_zero': False})
        outputs += [{'type': 'branch_current', 'component': fid}, {'type': 'bus_voltage', 'bus': bus_id(f['bus'])}]
    if case.get('record_currents'):
        outputs += [{'type': 'branch_current', 'component': c} for c in currents]
    return {'format': 'openkeito-model', 'version': FORMAT_VERSION, 'name': name,
            'description': 'IEEE 39-bus system with a WECC generic PV plant at bus 30, converted from the emthub CIM '
                           '(IEEE39.xml, IEEE39_ic.xml) by cim_to_openkeito.py',
            'simulation': {'time_step_s': dt, 'end_time_s': float(case['end_time_s']), 'method': method,
                           'output_interval_steps': int(case['output_interval_steps']), 'initialization': 'steady_state'},
            'components': comps, 'controls': controls, 'outputs': outputs}


# ---------------------------------------------------------------- sequence-impedance checks (frequency scan)

def _scan_model(name, comps):
    return {'format': 'openkeito-model', 'version': FORMAT_VERSION, 'name': f'scan_{name}',
            'simulation': {'time_step_s': 1e-5, 'end_time_s': 1e-4},
            'components': comps, 'frequency_scan': {'frequencies_hz': [FREQ], 'buses': ['P']}}


def _terminate(bus, r):
    return [{'type': 'resistor', 'name': f'R{bus}{p}', 'nodes': [f'{bus}.{p}', 'gnd'], 'resistance_ohm': r} for p in 'abc']


def _at(c, mapping):
    c = json.loads(json.dumps(c))
    if 'buses' in c:
        c['buses'] = [mapping.get(b, b) for b in c['buses']]
    if 'bus' in c:
        c['bus'] = mapping.get(c['bus'], c['bus'])
    return c


def sequence_tests(net, lines=('15_16_1', '16_19_1'), xfmrs=('19_33_0_1', '19_20_0_1')):
    """Small models for ``openkeito scan`` at 60 Hz (bus P) and the sequence impedances [ohm] expected from the CIM:
    lines with the far end shorted or open (exact distributed values Zc·tanh(gl), Zc·coth(gl)); transformers from the
    high-voltage side with the other side shorted ((r + jx)·tap², and for a delta secondary also open in the zero
    sequence); machines (negative sequence Ra + j(Xd'' + Xq'')/2, zero sequence Ra + jXl) and the swing source."""
    short, open_ = 1e-6, 1e9
    out = []
    by = {ln['name']: ln for ln in net['lines']}
    for nm in lines:
        ln = by[nm]
        c = line(ln, 1e-5)[0]
        c = _at(c, {c['buses'][0]: 'P', c['buses'][1]: 'Q'})
        exp = {}
        for seq, (r, x, b) in (('Z1', (ln['r'], ln['x'], ln['bch'])), ('Z0', (ln['r0'], ln['x0'], ln['b0ch']))):
            zc, gl = cmath.sqrt(complex(r, x) / complex(0, b)), cmath.sqrt(complex(r, x) * complex(0, b))
            exp[seq] = (zc * cmath.tanh(gl), zc / cmath.tanh(gl))
        out.append((f'line_{nm}_short', _scan_model(nm, [c] + _terminate('Q', short)), {k: v[0] for k, v in exp.items()}))
        out.append((f'line_{nm}_open', _scan_model(nm, [c]), {k: v[1] for k, v in exp.items()}))
    bt = {t['name']: t for t in net['xfmrs']}
    for nm in xfmrs:
        t = bt[nm]
        c = transformer(t)
        c = _at(c, {c['buses'][0]: 'P', c['buses'][1]: 'Q'})
        z = complex(t['r'], t['x']) * t['tap'] ** 2
        out.append((f'xfmr_{nm}_short', _scan_model(nm, [c] + _terminate('Q', short)), {'Z1': z, 'Z2': z, 'Z0': z}))
        if c['vector_group'].endswith('d1'):
            out.append((f'xfmr_{nm}_open', _scan_model(nm, [c] + _terminate('Q', open_)), {'Z0': z}))
    for g in net['gens']:
        zb = (g['ratedU'] / 1e3) ** 2 / (g['ratedS'] / 1e6)
        if g['bus'] == SWING_BUS:
            c = _at(swing_source(g, net), {bus_id(g['bus']): 'P'})
            c['line_voltage_rms_v'] = 0.0
            z = complex(c['impedance']['r1_ohm'], W * c['impedance']['l1_h'])
            out.append((f"swing_{g['name']}", _scan_model(g['name'], [c]), {'Z1': z, 'Z2': z, 'Z0': z}))
            continue
        c = _at(machine(g, net), {bus_id(g['bus']): 'P'})
        z2 = complex(g['statorResistance'], (g['xDirectSubtrans'] + xq2(g)) / 2) * zb
        z0 = complex(g['statorResistance'], g['statorLeakageReactance']) * zb
        out.append((f"gen_{g['name']}", _scan_model(g['name'], [c]), {'Z2': z2, 'Z0': z0}))
    return out


def write_models(cim, ic, out, cases=None):
    net = network(load_cim(cim, ic))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for name in cases or CASES:
        model = build_model(net, CASES[name], f'ieee39_{name}')
        (out / f'{name}.json').write_text(json.dumps(model, indent=1) + '\n', encoding='utf-8', newline='\n')
    scan = out / 'scan'
    scan.mkdir(exist_ok=True)
    expected = {}
    for name, model, exp in sequence_tests(net):
        (scan / f'{name}.json').write_text(json.dumps(model, indent=1) + '\n', encoding='utf-8', newline='\n')
        expected[name] = {k: [v.real, v.imag] for k, v in exp.items()}
    (scan / 'expected_ohm.json').write_text(json.dumps(expected, indent=1) + '\n', encoding='utf-8', newline='\n')
    return net


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('cim')
    ap.add_argument('ic')
    ap.add_argument('out')
    ap.add_argument('--case', action='append', choices=sorted(CASES))
    a = ap.parse_args()
    write_models(a.cim, a.ic, a.out, a.case)


if __name__ == '__main__':
    main()
