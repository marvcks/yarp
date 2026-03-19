#!/usr/bin/env python
"""
Run autoCG in parallel on a list of atom-mapped reaction SMILES.

Input:  JSON file (flat list of mapped SMILES strings) from enumerate_rxn.py
Output: JSON file with conformer results.

Usage:
    python run_autocg.py -i results.json -o autocg_results.json -j 4
"""

import os
import sys
import json
import tempfile
import subprocess
import multiprocessing
from pathlib import Path
from typing import List, Tuple, Optional
from rdkit import Chem

_script_dir = Path(__file__).resolve().parent  # gluco_example/
_AUTOCG_PATH = Path(os.environ.get('AUTOCG_PATH', str(_script_dir.parent.parent / 'autoCG')))
_YARP_BASE = Path(os.environ.get('YARP_BASE', str(_script_dir.parent.parent)))


def _strip_atom_mapping(smiles: str) -> str:
    """Remove atom mapping numbers from a reaction SMILES."""
    parts = smiles.split('>>')
    if len(parts) != 2:
        return smiles
    result = []
    for part in parts:
        frags = part.split('.')
        stripped = []
        for frag in frags:
            try:
                mol = Chem.MolFromSmiles(frag)
                if mol is None:
                    stripped.append(frag)
                    continue
                for atom in mol.GetAtoms():
                    atom.SetAtomMapNum(0)
                stripped.append(Chem.MolToSmiles(mol))
            except Exception:
                stripped.append(frag)
        result.append('.'.join(stripped))
    return '>>'.join(result)


def _parse_good_conformers(guess_log: Path) -> List[int]:
    """Parse guess.log to get good conformer indices."""
    try:
        lines = guess_log.read_text().splitlines()
        for line in reversed(lines):
            if 'are expected to be good conformers' in line:
                idx_part = line.split('idx:')[1].split('are')[0].strip()
                import re
                return [int(x) for x in re.findall(r'\d+', idx_part)]
        return []
    except Exception:
        return []


def _run_one(args):
    """Worker function for multiprocessing. Returns list of result dicts."""
    i, mapped_smiles, output_dir, chg, mult, use_mapped, se, cs, timeout = args

    rxn_dir = Path(output_dir) / f"rxn_{i}"
    rxn_dir.mkdir(parents=True, exist_ok=True)
    abs_rxn_dir = rxn_dir.resolve()
    working_dir = Path(tempfile.mkdtemp(prefix=f"autocg_rxn{i}_"))
    working_dir.mkdir(parents=True, exist_ok=True)

    smiles_for_autocg = mapped_smiles if use_mapped else _strip_atom_mapping(mapped_smiles)

    cmd = [
        sys.executable,
        str(_AUTOCG_PATH / 'generate.py'),
        smiles_for_autocg,
        '-sd', str(abs_rxn_dir),
        '-wd', str(working_dir),
        '--chg', str(chg),
        '--mult', str(mult),
    ]
    if se:
        cmd += ['-se', str(se)]
    if cs:
        cmd += ['-cs', str(cs)]

    env = os.environ.copy()
    env['PYTHONPATH'] = str(_YARP_BASE) + ':' + env.get('PYTHONPATH', '')
    env.setdefault('CALCULATOR', 'orca')

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        if result.returncode != 0 and result.stderr:
            print(f"  [rxn_{i}] autoCG stderr: {result.stderr[-300:]}", file=sys.stderr)
    except subprocess.TimeoutExpired:
        print(f"  [rxn_{i}] autoCG timed out", file=sys.stderr)
        return []
    except Exception as e:
        print(f"  [rxn_{i}] autoCG exception: {e}", file=sys.stderr)
        return []

    guess_log = abs_rxn_dir / 'guess.log'
    good_idxs = _parse_good_conformers(guess_log) if guess_log.exists() else []

    results = []
    for idx in good_idxs:
        r_xyz = abs_rxn_dir / str(idx) / 'R.xyz'
        p_xyz = abs_rxn_dir / str(idx) / 'P.xyz'
        if r_xyz.exists() and p_xyz.exists():
            results.append({
                'smiles': mapped_smiles,
                'rxn_dir': str(abs_rxn_dir),
                'conformer': idx,
                'r_xyz': str(r_xyz),
                'p_xyz': str(p_xyz),
            })

    return results


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run autoCG in parallel on mapped reaction SMILES')
    parser.add_argument('-i', '--input', required=True, help='Input JSON file (list of mapped SMILES)')
    parser.add_argument('-o', '--output', default=None, help='Output JSON file (default: autocg_results.json next to input)')
    parser.add_argument('-d', '--output-dir', default=None, dest='output_dir',
                        help='Directory for autoCG output subdirs (default: same dir as input JSON)')
    parser.add_argument('-j', '--jobs', type=int, default=1, help='Parallel processes (default: 1)')
    parser.add_argument('--chg', type=int, default=0, help='Total charge (default: 0)')
    parser.add_argument('--mult', type=int, default=1, help='Spin multiplicity (default: 1)')
    parser.add_argument('--no-mapping', action='store_true', help='Strip atom mapping before passing to autoCG')
    parser.add_argument('-se', '--stereo-enumerate', type=int, default=0, dest='se',
                        help='Whether to consider stereoisomers (default: 0)')
    parser.add_argument('-cs', '--check-stereo', type=int, default=0, dest='cs',
                        help='Screening criteria for removing same pseudo-TS (default: 0)')
    parser.add_argument('--timeout', type=int, default=3600, help='Per-reaction timeout in seconds (default: 3600)')

    args = parser.parse_args()

    input_path = Path(args.input)
    with open(input_path) as f:
        smiles_list = json.load(f)

    output_dir = Path(args.output_dir) if args.output_dir else input_path.parent / 'autocg_runs'
    output_dir.mkdir(parents=True, exist_ok=True)

    json_out = args.output if args.output else str(input_path.parent / 'autocg_results.json')

    use_mapped = not args.no_mapping

    tasks = [
        (i, smi, str(output_dir), args.chg, args.mult, use_mapped, args.se, args.cs, args.timeout)
        for i, smi in enumerate(smiles_list)
    ]

    if args.jobs == 1:
        all_results = []
        for task in tasks:
            all_results.extend(_run_one(task))
    else:
        workers = args.jobs if args.jobs > 0 else multiprocessing.cpu_count()
        with multiprocessing.Pool(workers) as pool:
            nested = pool.map(_run_one, tasks)
        all_results = [r for sublist in nested for r in sublist]

    with open(json_out, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(json_out)
