#!/usr/bin/env python
"""
Enumerate reactions from a reactant XYZ file, output atom-mapped reaction SMILES list (JSON).

Usage:
    python enumerate_rxn.py -i input.xyz -o output_dir
    python enumerate_rxn.py -i input.xyz -o output_dir -n 2 --max 10
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
from typing import List, Optional, Tuple

# Add paths - configurable via environment variables
_script_dir = Path(__file__).resolve().parent  # gluco_example/
_yarp_dir = Path(os.environ.get('YARP_BASE', str(_script_dir.parent.parent)))  # two levels up by default
_erm_dir = Path(os.environ.get('ERM_PATH', str(_script_dir.parent / 'Enumerate-Reaction-Mapping-main')))
sys.path.insert(0, str(_yarp_dir))
sys.path.insert(0, str(_erm_dir))

import yarp as yp
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from yarp.find_lewis import return_formals
from enumerate_reaction_mapping import generate_atom_mapped_reactions


def return_smi(elements, bond_mat, fc) -> Tuple[str, str]:
    """Generate atom-mapped and plain SMILES from elements/bond_mat/fc."""
    try:
        mol = Chem.RWMol()
        for i, element in enumerate(elements):
            atom = Chem.Atom(element.capitalize())
            atom.SetAtomMapNum(i + 1)
            atom.SetFormalCharge(int(fc[i]))
            atom.SetNumRadicalElectrons(int(bond_mat[i, i] % 2))
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
            mol.AddAtom(atom)
        num_atoms = len(elements)
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                order = int(round(bond_mat[i, j]))
                if order == 1:
                    bond_type = Chem.rdchem.BondType.SINGLE
                elif order == 2:
                    bond_type = Chem.rdchem.BondType.DOUBLE
                elif order == 3:
                    bond_type = Chem.rdchem.BondType.TRIPLE
                elif order == 4:
                    bond_type = Chem.rdchem.BondType.QUADRUPLE
                else:
                    continue
                mol.AddBond(i, j, bond_type)
        Chem.SanitizeMol(mol)
        mapped_smiles = Chem.MolToSmiles(mol)
        plain_mol = Chem.RWMol(mol)
        for atom in plain_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        plain_smiles = Chem.MolToSmiles(plain_mol)
        return mapped_smiles, plain_smiles
    except Exception:
        return "", ""


def get_bond_changes(reactant_bond_mat, product_bond_mat) -> Tuple[List, List]:
    """Compare bond matrices to identify bond formation and breaking."""
    n = len(reactant_bond_mat)
    formed = []
    broken = []
    for i in range(n):
        for j in range(i + 1, n):
            r_bond = reactant_bond_mat[i, j]
            p_bond = product_bond_mat[i, j]
            if p_bond > r_bond:
                formed.append((i, j))
            elif r_bond > p_bond:
                broken.append((i, j))
    return formed, broken


def enumerate_rxn(
    r_xyz: str,
    n_break: int = 1,
    output_dir: Optional[str] = None,
    chg: int = 0,
    max_reactions: Optional[int] = None,
    verbose: bool = False,
) -> List[str]:
    """
    Enumerate reactions from a reactant XYZ file.
    Returns a flat list of atom-mapped reaction SMILES strings.
    """
    mol = yp.yarpecule(r_xyz)
    if verbose:
        print(f"Input: {r_xyz}", file=sys.stderr)
        print(f"Atoms: {len(mol.elements)}", file=sys.stderr)

    if output_dir is None:
        output_dir = Path(r_xyz).parent / "reactions"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    hashes = set([mol.hash])
    broken = list(set(yp.break_bonds([mol], n=n_break)))
    if verbose:
        print(f"Intermediates after breaking {n_break} bond(s): {len(broken)}", file=sys.stderr)

    all_products = []
    for intermediate in broken:
        products = list(set(yp.form_bonds([intermediate], hashes=hashes, inter=True)))
        all_products.extend(products)
    all_products = list(set(all_products))
    if verbose:
        print(f"Products: {len(all_products)}", file=sys.stderr)

    reactant_bond_mat = np.array(mol.bond_mats[0])
    reactant_elements = [e.upper() for e in mol.elements]
    react_formal_charge = return_formals(reactant_bond_mat, reactant_elements)

    reaction_pairs = []
    for idx, product in enumerate(all_products):
        product_bond_mat = np.array(product.bond_mats[0])
        product_elements = [e.upper() for e in product.elements]
        prod_formal_charge = return_formals(product_bond_mat, [e.lower() for e in product.elements])

        r_mapped, _ = return_smi(reactant_elements, reactant_bond_mat, react_formal_charge)
        p_mapped, _ = return_smi(product_elements, product_bond_mat, prod_formal_charge)

        if not r_mapped or not p_mapped:
            continue

        formed, broken_bonds = get_bond_changes(reactant_bond_mat, product_bond_mat)
        if len(formed) == 0 and len(broken_bonds) == 0:
            continue

        reaction_pairs.append((idx, r_mapped, p_mapped))

    if verbose:
        print(f"Valid reaction pairs: {len(reaction_pairs)}", file=sys.stderr)

    results = []
    processed_count = 0

    for idx, r_mapped, p_mapped in reaction_pairs:
        if max_reactions is not None and processed_count >= max_reactions:
            if verbose:
                print(f"Reached max_reactions limit ({max_reactions})", file=sys.stderr)
            break

        mapped_rxn = f"{r_mapped}>>{p_mapped}"
        try:
            mapped_smiles_list = generate_atom_mapped_reactions(rxn_smiles=mapped_rxn, debug=False, jobs=1)
        except Exception as e:
            print(f"  [{idx}] atom mapping failed: {e}", file=sys.stderr)
            continue

        if not mapped_smiles_list:
            continue

        for mapped_smiles in mapped_smiles_list:
            results.append(mapped_smiles)
            processed_count += 1

    if verbose:
        print(f"Total mapped reactions: {len(results)}", file=sys.stderr)

    return results


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Enumerate reactions from a reactant XYZ file')
    parser.add_argument('-i', '--input', required=True, help='Reactant XYZ file path')
    parser.add_argument('-o', '--output', default=None, help='Output directory')
    parser.add_argument('-n', '--break-n', type=int, default=1, dest='break_n',
                        help='Number of bonds to break (default: 1)')
    parser.add_argument('--max', dest='max_reactions', type=int, default=None,
                        help='Max number of reactions to enumerate')
    parser.add_argument('--chg', type=int, default=0, help='Total charge (default: 0)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output to stderr')
    parser.add_argument('--json', dest='json_output', default=None,
                        help='Output JSON file path (default: output_dir/results.json)')

    args = parser.parse_args()

    results = enumerate_rxn(
        r_xyz=args.input,
        n_break=args.break_n,
        output_dir=args.output,
        chg=args.chg,
        max_reactions=args.max_reactions,
        verbose=args.verbose,
    )

    output_dir = Path(args.output) if args.output else Path(args.input).parent / 'reactions'
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.json_output if args.json_output else str(output_dir / 'results.json')
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(json_path)
