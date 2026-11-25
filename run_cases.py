# SPDX-License-Identifier: MIT OR Apache-2.0
# Copyright (c) 2026 The OpenKeito contributors
"""Run the eight IEEE 39-bus cases and the sequence-impedance scans with OpenKeito.

    python run_cases.py [--openkeito PATH] [--timing] [CASE ...]

For every case ``models/<case>.json`` this writes ``runs/<case>/`` with ``version.txt`` (``openkeito
--version``), ``validate.json`` and ``run/`` (``result.csv``, ``phasors.csv``, ``events.csv``,
``summary.json``).  The scan models ``models/scan/*.json`` go to ``runs/scan/<name>/`` and are compared with
the impedances expected from the CIM (``models/scan/expected_ohm.json``) by ``make_figures.py``.

``--timing`` instead measures the computing time: ``nofault``, ``slg16`` and ``slg16_20s`` are run three times
each with ``--threads 1`` on one logical core (``start /affinity 1`` on Windows, ``sched_setaffinity`` on
Linux), and the minimum of ``elapsed_s`` of ``summary.json`` is written to ``runs/timing.json``.

Only the Python standard library is needed.  The OpenKeito command is taken from ``--openkeito``, else from
the environment variable ``OPENKEITO``, else ``openkeito`` on the PATH.
"""
import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = ['nofault', 'nofault_qreg0', 'slg16', 'tpf02', 'slg16_20s', 'nofault_conditional_pi',
         'slg16_20s_mag_transformer', 'slg16_20s_mag_shunt']
TIMING_CASES = ['nofault', 'slg16', 'slg16_20s']


def openkeito(exe, *args, stdout=None):
    r = subprocess.run([exe, *map(str, args)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding='utf-8')
    if stdout is not None:
        Path(stdout).write_text(r.stdout, encoding='utf-8', newline='\n')
    return r


def run_case(exe, case, runs):
    d = runs / case
    d.mkdir(parents=True, exist_ok=True)
    model = HERE / 'models' / f'{case}.json'
    openkeito(exe, '--version', stdout=d / 'version.txt')
    v = openkeito(exe, 'validate', model, stdout=d / 'validate.json')
    r = openkeito(exe, 'run', model, '--out', d / 'run', '--force', stdout=d / 'run.json')
    if r.returncode != 0:
        sys.exit(f'{case}: openkeito run failed (exit {r.returncode})\n{r.stdout}\n{r.stderr}')
    s = json.loads((d / 'run' / 'summary.json').read_text(encoding='utf-8'))
    val = json.loads(v.stdout)
    codes = sorted({w['code'] for w in val['warnings'] + s['warnings']})
    print(f"{case}: valid {val['valid']}, {s['status']}, {s['elapsed_s']:.1f} s, warnings {codes}")


def run_scans(exe, runs):
    """openkeito scan of every model in models/scan into runs/scan/<name>/ (make_figures.py compares them with the
    impedances expected from the CIM)."""
    src = HERE / 'models' / 'scan'
    expected = json.loads((src / 'expected_ohm.json').read_text(encoding='utf-8'))
    for name in expected:
        r = openkeito(exe, 'scan', src / f'{name}.json', '--out', runs / 'scan' / name, '--force')
        if r.returncode != 0:
            sys.exit(f'scan {name} failed (exit {r.returncode})\n{r.stdout}\n{r.stderr}')
    print(f'scan: {len(expected)} models')


def machine():
    if sys.platform == 'win32':
        ps = ('(Get-CimInstance Win32_Processor).Name; (Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory; '
              '[Environment]::OSVersion.VersionString')
        out = subprocess.run(['powershell', '-NoProfile', '-Command', ps], capture_output=True, text=True).stdout.splitlines()
        if len(out) >= 3:
            return out[0].strip(), int(out[1]), out[2].strip()
    cpu = platform.processor()
    try:
        for line in Path('/proc/cpuinfo').read_text().splitlines():
            if line.startswith('model name'):
                cpu = line.split(':', 1)[1].strip()
                break
    except OSError:
        pass
    mem = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') if hasattr(os, 'sysconf') else None
    return cpu, mem, platform.platform()


def run_timing(exe, runs, repeat=3):
    if sys.platform == 'win32':
        affinity = 'one logical core (start /affinity 1)'
        exe_path = shutil.which(exe) or exe
    elif hasattr(os, 'sched_setaffinity'):
        os.sched_setaffinity(0, {sorted(os.sched_getaffinity(0))[0]})
        affinity = 'one logical core (sched_setaffinity)'
    else:
        affinity = 'not pinned'
    cpu, mem, osname = machine()
    out = {'cpu': cpu, 'memory_bytes': mem, 'os': osname, 'affinity': affinity, 'threads': 1, 'repeat': repeat,
           'statistic': 'minimum of elapsed_s over the repeats', 'cases': {}}
    for case in TIMING_CASES:
        d = runs / 'timing' / case
        model = HERE / 'models' / f'{case}.json'
        times = []
        for _ in range(repeat):
            args = ['run', str(model), '--out', str(d), '--force', '--threads', '1']
            if sys.platform == 'win32':
                cmd = 'start "" /wait /b /affinity 1 ' + subprocess.list2cmdline([exe_path, *args])
                subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, check=True)
            else:
                subprocess.run([exe, *args], stdout=subprocess.DEVNULL, check=True)
            s = json.loads((d / 'summary.json').read_text(encoding='utf-8'))
            times.append(s['elapsed_s'])
        steps = s['solver']['steps']
        out['cases'][case] = {'elapsed_s': min(times), 'all_elapsed_s': times, 'steps': steps,
                              'us_per_step': 1e6 * min(times) / steps}
        print(case, out['cases'][case])
    (runs / 'timing.json').write_text(json.dumps(out, indent=1) + '\n', encoding='utf-8', newline='\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('cases', nargs='*', help=f"cases to run (default: all eight and the scans): {', '.join(CASES)}")
    ap.add_argument('--openkeito', default=os.environ.get('OPENKEITO', 'openkeito'), help='the openkeito command')
    ap.add_argument('--runs', default=str(HERE / 'runs'), help='output folder (default: runs)')
    ap.add_argument('--timing', action='store_true', help='measure the computing time instead')
    a = ap.parse_args()
    bad = [c for c in a.cases if c not in CASES]
    if bad:
        ap.error(f'unknown case {bad}')
    runs = Path(a.runs)
    runs.mkdir(parents=True, exist_ok=True)
    if a.timing:
        return run_timing(a.openkeito, runs)
    for case in a.cases or CASES:
        run_case(a.openkeito, case, runs)
    if not a.cases:
        run_scans(a.openkeito, runs)


if __name__ == '__main__':
    main()
