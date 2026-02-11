#!/usr/bin/env python3
import json
import shutil
import subprocess
from pathlib import Path

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / 'outputs' / 'phase_constraints'
INPUT_DIR = REPO_ROOT / 'inputs' / 'phase_constraints'
DEFAULT_JSON = REPO_ROOT / 'outputs' / 'thermochimica_out.json'

TEMP_SWEEP_CASES = [
    ('unconstrained', INPUT_DIR / 'co_temp_sweep_unconstrained.ti'),
    ('constrained 50/50', INPUT_DIR / 'co_temp_sweep_constrained_50_50.ti'),
    ('constrained 70/30', INPUT_DIR / 'co_temp_sweep_constrained_70_30.ti'),
]

COMP_SWEEP_CASES = [
    ('unconstrained', INPUT_DIR / 'co_comp_sweep_unconstrained.ti'),
    ('constrained 50/50', INPUT_DIR / 'co_comp_sweep_constrained_50_50.ti'),
    ('constrained 70/30', INPUT_DIR / 'co_comp_sweep_constrained_70_30.ti'),
]

PHASES = ['gas_ideal', 'C_Graphite(s)']


def run_input_script(input_path: Path, json_out: Path):
    subprocess.run([str(REPO_ROOT / 'bin' / 'InputScriptMode'), str(input_path)], check=True)
    shutil.copy2(DEFAULT_JSON, json_out)


def run_calc_list(input_path: Path, json_out: Path):
    subprocess.run([str(REPO_ROOT / 'bin' / 'RunCalculationList'), str(input_path)], check=True)
    shutil.copy2(DEFAULT_JSON, json_out)


def load_series(json_file: Path, phase_name: str):
    with open(json_file, 'r') as f:
        data = json.load(f)
    keys = sorted(data.keys(), key=lambda k: int(k))
    temps = []
    ratios = []
    fracs = []
    for k in keys:
        entry = data[k]
        if not entry or 'elements' not in entry:
            continue
        total_elems = sum(e['moles'] for e in entry['elements'].values())
        if total_elems <= 0:
            total_elems = 1.0

        phase = None
        if phase_name in entry.get('solution phases', {}):
            phase = entry['solution phases'][phase_name]
        elif phase_name in entry.get('pure condensed phases', {}):
            phase = entry['pure condensed phases'][phase_name]

        phase_total = 0.0
        if phase is not None:
            phase_total = sum(e['moles of element in phase'] for e in phase['elements'].values())

        fracs.append(phase_total / total_elems)
        temps.append(entry['temperature'])
        if 'C' in entry['elements'] and 'O' in entry['elements']:
            c_moles = entry['elements']['C']['moles']
            o_moles = entry['elements']['O']['moles']
            ratio = o_moles / c_moles if c_moles > 0 else 0.0
        else:
            ratio = 0.0
        ratios.append(ratio)

    return temps, ratios, fracs


def find_missing_entries(json_file: Path):
    with open(json_file, 'r') as f:
        data = json.load(f)
    keys = sorted(data.keys(), key=lambda k: int(k))
    return [k for k in keys if not data[k] or 'elements' not in data[k]]


def write_csv(results, xkey, out_path: Path):
    header = ['x', 'phase', 'case', 'phase_fraction']
    lines = [','.join(header)]
    for label, json_file in results:
        for phase in PHASES:
            temps, ratios, fracs = load_series(json_file, phase)
            xs = temps if xkey == 'temperature' else ratios
            for x, frac in zip(xs, fracs):
                lines.append(f'{x},{phase},{label},{frac}')
    out_path.write_text('\n'.join(lines) + '\n')


def plot_temp_sweep(results, output_path: Path):
    if not HAS_MPL:
        return
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8.0), sharex=True)
    for label, json_file in results:
        for idx, phase in enumerate(PHASES):
            temps, _, fracs = load_series(json_file, phase)
            axes[idx].plot(temps, fracs, label=label)
            axes[idx].set_ylabel(f'{phase} phase fraction')
            axes[idx].set_ylim(0.0, 1.0)
            axes[idx].grid(True, alpha=0.3)
    axes[-1].set_xlabel('Temperature (K)')
    axes[0].legend(loc='best')
    fig.suptitle('CO system: phase fraction vs temperature')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_comp_sweep(results, output_path: Path):
    if not HAS_MPL:
        return
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8.0), sharex=True)
    for label, json_file in results:
        for idx, phase in enumerate(PHASES):
            _, ratios, fracs = load_series(json_file, phase)
            axes[idx].plot(ratios, fracs, label=label)
            axes[idx].set_ylabel(f'{phase} phase fraction')
            axes[idx].set_ylim(0.0, 1.0)
            axes[idx].grid(True, alpha=0.3)
    axes[-1].set_xlabel('O/C (mole ratio)')
    axes[0].legend(loc='best')
    fig.suptitle('CO system (1000 K): phase fraction vs O/C')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    temp_results = []
    for tag, input_path in TEMP_SWEEP_CASES:
        json_out = OUTPUT_DIR / f'co_temp_{tag.replace(" ", "_").replace("/", "-")}.json'
        run_input_script(input_path, json_out)
        temp_results.append((tag, json_out))

    comp_results = []
    for tag, input_path in COMP_SWEEP_CASES:
        json_out = OUTPUT_DIR / f'co_comp_{tag.replace(" ", "_").replace("/", "-")}.json'
        run_calc_list(input_path, json_out)
        comp_results.append((tag, json_out))

    for _, json_file in temp_results + comp_results:
        missing = find_missing_entries(json_file)
        if missing:
            print(f'Warning: {json_file.name} missing results for steps: {", ".join(missing)}')

    write_csv(temp_results, 'temperature', OUTPUT_DIR / 'co_temp_sweep_phase_fractions.csv')
    write_csv(comp_results, 'ratio', OUTPUT_DIR / 'co_comp_sweep_phase_fractions.csv')

    plot_temp_sweep(temp_results, OUTPUT_DIR / 'co_temp_sweep_phase_fractions.png')
    plot_comp_sweep(comp_results, OUTPUT_DIR / 'co_comp_sweep_phase_fractions.png')

    if HAS_MPL:
        print('Wrote plots to:', OUTPUT_DIR)
    else:
        print('matplotlib not available; wrote CSVs to:', OUTPUT_DIR)
        print('Install matplotlib and rerun to generate PNG plots.')


if __name__ == '__main__':
    main()
