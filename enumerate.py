#!/usr/bin/env python
"""
Enumerate bNfN reactions using YARP from a SMILES input.

Usage:
    python enumerate.py -s <SMILES> -m b2f2 -o output.json
    python enumerate.py -s <SMILES> -m b2f2 -o output.json --break-higher-order
"""
import re
import json
import argparse

import numpy as np

from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import yarp as yp
from yarp.find_lewis import return_formals


# Default HMF target SMILES (canonicalized)
DEFAULT_TARGET_SMILES = "O=Cc1ccc(CO)o1.O.O.O"


def canonicalize_smiles(smiles: str) -> str | None:
    """Return canonical SMILES with stereochemistry removed."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)


def canonicalize_product_smiles(reaction_smiles: str) -> str | None:
    """Extract and canonicalize the product side of a reaction SMILES."""
    parts = reaction_smiles.split(">>", 1)
    if len(parts) != 2:
        return None
    return canonicalize_smiles(parts[1])


def check_target_in_reactions(reactions: list, target_smiles: str) -> bool:
    """Check if target SMILES exists as a product in any reaction."""
    target_canon = canonicalize_smiles(target_smiles)
    if target_canon is None:
        return False
    for rxn in reactions:
        product_canon = canonicalize_product_smiles(rxn)
        if product_canon == target_canon:
            return True
    return False


def get_fragments(side):
    """Return list of fragment SMILES from one side of a reaction SMILES."""
    return side.split('.')


def get_product_fragments(rxn_smi):
    """Return list of fragment SMILES from the product side."""
    parts = rxn_smi.split('>>')
    if len(parts) != 2:
        return []
    return parts[1].split('.')


def is_water(fragment):
    """Match water in any common SMILES representation: [H]O[H], O, [OH2], etc."""
    f = fragment.strip()
    patterns = [
        r'\[H(:\d+)?\]O\[H(:\d+)?\]',
        r'\[H(:\d+)?\]\[O(:\d+)?\]\[H(:\d+)?\]',
        r'O',
        r'\[OH2(:\d+)?\]',
        r'\[O@@H2\]',
        r'\[O@H2\]',
    ]
    return any(re.fullmatch(p, f) for p in patterns)


def filter_reactions(reaction_smiles_list):
    """Keep reactions where the only released small molecules are water.

    A reaction is kept if all product fragments except at most one are water.
    Any non-water small molecule in the products causes the reaction to be dropped.
    """
    filtered = []
    for smi in reaction_smiles_list:
        parts = smi.split('>>')
        if len(parts) != 2:
            continue
        product_frags = get_fragments(parts[1])
        non_water_frags = [f for f in product_frags if not is_water(f)]
        if len(non_water_frags) <= 1:
            filtered.append(smi)
    return filtered


def bond_mat_to_smiles(elements, bond_mat, fc):
    """Generate SMILES from elements, bond_mat, and formal charges."""
    try:
        mol = Chem.RWMol()
        for i, element in enumerate(elements):
            atom = Chem.Atom(element.capitalize())
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
                else:
                    continue
                mol.AddBond(i, j, bond_type)

        Chem.SanitizeMol(mol)
        smiles = Chem.MolToSmiles(mol)
        return smiles
    except Exception:
        return None


def enumerate_reactions(smiles, n_break, output_file=None, break_higher_order=False):
    """
    Enumerate bNfN reactions from a SMILES string.

    Parameters
    ----------
    smiles : str
        Input SMILES string (can be reaction SMILES with >>)
    n_break : int
        Number of bonds to break (and form), e.g. 1, 2, or 3
    output_file : str, optional
        Output JSON file path
    break_higher_order : bool, optional
        Whether to break higher-order bonds (double/triple) (default: False)

    Returns
    -------
    reaction_smiles_list : list
        List of reaction SMILES strings
    """
    # Handle reaction SMILES format: extract reactant side
    if '>>' in smiles:
        parts = smiles.split('>>')
        reactant_smiles = parts[0]
    else:
        reactant_smiles = smiles
        
    reactant_smiles = canonicalize_smiles(reactant_smiles)
    mol = yp.yarpecule(reactant_smiles)
    hashes = set([mol.hash])

    broken = list(set(yp.break_bonds([mol], n=n_break, break_higher_order=break_higher_order)))
    all_products = []
    for i, intermediate in enumerate(broken):
        products = list(
                    set(
                        yp.form_bonds(
                            [intermediate], 
                            hashes=hashes, 
                            inter=True,
                            intra=True,
                            def_only=True,
                            )
                        )
                    )
        all_products.extend(products)

    all_products = list(set(all_products))
    unique_products = [p for p in all_products if p.hash != mol.hash]

    # First pass: generate all reaction SMILES and save XYZ files
    all_reactions = []
    for idx, prod in enumerate(unique_products):
        prod_bond_mat = np.array(prod.bond_mats[0])
        prod_elements = [e.lower() for e in prod.elements]
        prod_fc = return_formals(prod_bond_mat, prod_elements)
        prod_smiles = bond_mat_to_smiles(prod_elements, prod_bond_mat, prod_fc)

        if prod_smiles:
            canon_smiles = canonicalize_smiles(prod_smiles)
            if canon_smiles:
                rxn_smiles = f"{reactant_smiles}>>{canon_smiles}"
                all_reactions.append({'rxn_smiles': rxn_smiles, 'prod': prod, 'idx': idx})

    # Filter reactions based on water criteria
    reaction_smiles_list = filter_reactions([r['rxn_smiles'] for r in all_reactions])

    if output_file:
        with open(output_file, 'w') as f:
            json.dump(reaction_smiles_list, f, indent=2)
        print(f"\nSaved to: {output_file}")

    return reaction_smiles_list

MODE_TO_N = {'b1f1': 1, 'b2f2': 2, 'b3f3': 3}

def main():
    parser = argparse.ArgumentParser(
        description='Enumerate bNfN reactions from a SMILES input using YARP.'
    )
    parser.add_argument('-s', '--smiles', required=True,
                        help='Input SMILES string')
    parser.add_argument('-m', '--mode', choices=['b1f1', 'b2f2', 'b3f3'], default='b2f2',
                        help='Reaction mode: b1f1, b2f2, or b3f3 (default: b2f2)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output JSON file path')
    parser.add_argument('--break-higher-order', action='store_true',
                        help='Break higher-order bonds (double/triple bonds)')
    parser.add_argument('-t', '--target', default=None,
                        help=f'Target product SMILES to check (default: {DEFAULT_TARGET_SMILES})')
    args = parser.parse_args()

    n_break = MODE_TO_N[args.mode]
    reactions = enumerate_reactions(args.smiles, n_break, args.output, break_higher_order=args.break_higher_order)

    # Check for target product
    target = args.target if args.target else DEFAULT_TARGET_SMILES
    if check_target_in_reactions(reactions, target):
        print(f"Target product '{target}' exists in reactions.")
    else:
        print(f"Target product '{target}' NOT found in reactions.")


if __name__ == '__main__':
    main()
