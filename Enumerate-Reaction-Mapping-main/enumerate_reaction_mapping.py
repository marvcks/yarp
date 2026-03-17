import argparse
import math
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

import networkx as nx
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

EPS = 1e-6


def parse_reaction_smiles(rxn_smiles):
    if ">>" in rxn_smiles:
        parts = rxn_smiles.split(">>")
        if len(parts) != 2:
            raise ValueError("Reaction SMILES should contain a single '>>'.")
        reactants, products = parts
    else:
        parts = rxn_smiles.split(">")
        if len(parts) != 3:
            raise ValueError("Reaction SMILES should contain '>>' or three '>' parts.")
        reactants, products = parts[0], parts[2]

    if not reactants or not products:
        raise ValueError("Reaction SMILES must include reactants and products.")

    return reactants, products


def smiles_to_mol(smiles):
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
        mol = Chem.AddHs(mol)
    except Exception as e:
        raise ValueError(f"Could not kekulize molecule from SMILES: {smiles}") from e
    return mol


def mol_to_graph(mol):
    graph = nx.Graph()
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx(), symbol=atom.GetSymbol())
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        graph.add_edge(a1, a2, bond_order=bond.GetBondTypeAsDouble())
    return graph


def atom_symbol_counts(mol):
    return Counter(atom.GetSymbol() for atom in mol.GetAtoms())


def symmetry_ranks(mol):
    return Chem.CanonicalRankAtoms(mol, breakTies=False, includeAtomMaps=False)


def mapping_key(atom_mapping, reactant_ranks, product_ranks):
    counts = {}
    for react_idx, prod_idx in atom_mapping.items():
        r_rank = reactant_ranks[react_idx]
        p_rank = product_ranks[prod_idx]
        counts.setdefault(r_rank, {})
        counts[r_rank][p_rank] = counts[r_rank].get(p_rank, 0) + 1

    key_parts = []
    for r_rank in sorted(set(reactant_ranks)):
        p_counts = counts.get(r_rank, {})
        key_parts.append((r_rank, tuple(sorted(p_counts.items()))))
    return tuple(key_parts)


def _atom_mapnum_to_idx(mol):
    mapnum_to_idx = {}
    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if not map_num:
            continue
        if map_num in mapnum_to_idx:
            return None
        mapnum_to_idx[map_num] = atom.GetIdx()
    return mapnum_to_idx


def _atom_mapnum_to_symbol(mol):
    mapnum_to_symbol = {}
    duplicates = set()
    for atom in mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if not map_num:
            continue
        if map_num in mapnum_to_symbol:
            duplicates.add(map_num)
        mapnum_to_symbol[map_num] = atom.GetSymbol()
    return mapnum_to_symbol, duplicates


def _describe_seed_skip_reason(reactant_mol, product_mol):
    react_map, react_dups = _atom_mapnum_to_symbol(reactant_mol)
    prod_map, prod_dups = _atom_mapnum_to_symbol(product_mol)

    if not react_map and not prod_map:
        return "no atom-map numbers found"

    if react_dups or prod_dups:
        react_preview = sorted(react_dups)[:10]
        prod_preview = sorted(prod_dups)[:10]
        return f"duplicate atom-map numbers (reactant={react_preview}, product={prod_preview})"

    only_react = sorted(set(react_map) - set(prod_map))
    only_prod = sorted(set(prod_map) - set(react_map))
    if only_react or only_prod:
        react_preview = only_react[:10]
        prod_preview = only_prod[:10]
        return f"atom-map sets differ (only reactant={react_preview}, only product={prod_preview})"

    mismatches = []
    for map_num in sorted(react_map):
        if react_map[map_num] != prod_map[map_num]:
            mismatches.append((map_num, react_map[map_num], prod_map[map_num]))
    if mismatches:
        preview = mismatches[:10]
        suffix = "" if len(mismatches) <= 10 else f" (+{len(mismatches) - 10} more)"
        return f"atom element mismatch for map nums: {preview}{suffix}"

    return "seed mapping could not be constructed"


def _mapping_from_atom_maps(reactant_mol, product_mol):
    react_map = _atom_mapnum_to_idx(reactant_mol)
    prod_map = _atom_mapnum_to_idx(product_mol)
    if not react_map or not prod_map:
        return None
    if set(react_map) != set(prod_map):
        return None

    mapping = {}
    for map_num, react_idx in react_map.items():
        prod_idx = prod_map[map_num]
        if (
            reactant_mol.GetAtomWithIdx(react_idx).GetSymbol()
            != product_mol.GetAtomWithIdx(prod_idx).GetSymbol()
        ):
            return None
        mapping[react_idx] = prod_idx
    return mapping


def _complete_mapping_for_pruning(
    reactant_mol, product_mol, mapping, reactant_symbols, product_indices_by_symbol
):
    n_atoms = reactant_mol.GetNumAtoms()
    if n_atoms != product_mol.GetNumAtoms():
        return None

    used_product = [False] * n_atoms
    completed = {}
    for react_idx, prod_idx in mapping.items():
        if prod_idx < 0 or prod_idx >= n_atoms:
            return None
        if used_product[prod_idx]:
            return None
        completed[react_idx] = prod_idx
        used_product[prod_idx] = True

    if len(completed) == n_atoms:
        return completed

    # Prefer mapping unmapped hydrogens to hydrogens attached to already-mapped neighbors.
    for react_atom in reactant_mol.GetAtoms():
        react_idx = react_atom.GetIdx()
        if react_idx in completed:
            continue
        if react_atom.GetSymbol() != "H":
            continue
        neighbors = list(react_atom.GetNeighbors())
        if len(neighbors) != 1:
            continue
        react_nbr = neighbors[0].GetIdx()
        prod_nbr = completed.get(react_nbr)
        if prod_nbr is None:
            continue
        prod_atom = product_mol.GetAtomWithIdx(prod_nbr)
        for prod_neighbor in prod_atom.GetNeighbors():
            prod_idx = prod_neighbor.GetIdx()
            if prod_neighbor.GetSymbol() != "H":
                continue
            if used_product[prod_idx]:
                continue
            completed[react_idx] = prod_idx
            used_product[prod_idx] = True
            break

    for react_idx in range(n_atoms):
        if react_idx in completed:
            continue
        symbol = reactant_symbols[react_idx]
        candidates = product_indices_by_symbol.get(symbol)
        if not candidates:
            return None
        for prod_idx in candidates:
            if used_product[prod_idx]:
                continue
            completed[react_idx] = prod_idx
            used_product[prod_idx] = True
            break
        else:
            return None

    if len(completed) != n_atoms:
        return None
    return completed


def _seed_best_from_input_mapping(
    reactant_mol,
    product_mol,
    *,
    reactant_orders,
    product_orders,
    reactant_symbols,
    product_indices_by_symbol,
):
    mapping = _mapping_from_atom_maps(reactant_mol, product_mol)
    if mapping is None:
        return None, None
    completed = _complete_mapping_for_pruning(
        reactant_mol, product_mol, mapping, reactant_symbols, product_indices_by_symbol
    )
    if completed is None:
        return None, None

    used_product = [False] * product_mol.GetNumAtoms()
    tmp_mapping = {}
    add_count = 0.0
    break_count = 0.0
    for react_idx in range(reactant_mol.GetNumAtoms()):
        prod_idx = completed.get(react_idx)
        if prod_idx is None:
            return None, None
        if used_product[prod_idx]:
            return None, None
        da, db = _apply_assignment(
            tmp_mapping,
            used_product,
            react_idx,
            prod_idx,
            reactant_orders,
            product_orders,
        )
        add_count += da
        break_count += db

    best_total = add_count + break_count
    best_pairs = {(break_count, add_count)}
    return best_total, best_pairs


def bond_order_matrix(mol):
    n_atoms = mol.GetNumAtoms()
    matrix = [[0.0 for _ in range(n_atoms)] for _ in range(n_atoms)]
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        order = bond.GetBondTypeAsDouble()
        matrix[a1][a2] = order
        matrix[a2][a1] = order
    return matrix


def group_atoms_by_symbol(graph):
    groups = {}
    for node, data in graph.nodes(data=True):
        groups.setdefault(data["symbol"], []).append(node)
    return groups


def estimate_mapping_count(groups):
    total = 1
    for nodes in groups.values():
        total *= math.factorial(len(nodes))
    return total


def add_break_bounds_from_delta(min_delta, max_delta):
    if max_delta <= 0:
        min_add = 0.0
        max_add = 0.0
        min_break = -max_delta
        max_break = -min_delta
    elif min_delta >= 0:
        min_add = min_delta
        max_add = max_delta
        min_break = 0.0
        max_break = 0.0
    else:
        min_add = 0.0
        max_add = max_delta
        min_break = 0.0
        max_break = -min_delta
    return min_add, max_add, min_break, max_break


def format_count(value):
    if abs(value - round(value)) < EPS:
        return str(int(round(value)))
    return f"{value:.3f}".rstrip("0").rstrip(".")


def format_best_pairs(pairs):
    if not pairs:
        return "none"
    parts = []
    for break_count, add_count in sorted(pairs):
        parts.append(f"b{format_count(break_count)}f{format_count(add_count)}")
    return ",".join(parts)


def minmax_product_orders(candidates_a, candidates_b, product_orders, allow_same):
    min_p = None
    max_p = None
    for a in candidates_a:
        row = product_orders[a]
        for b in candidates_b:
            if not allow_same and a == b:
                continue
            val = row[b]
            if min_p is None or val < min_p:
                min_p = val
            if max_p is None or val > max_p:
                max_p = val
    return min_p, max_p


def compute_remaining_min_total(
    order,
    pos,
    mapping,
    used_product,
    reactant_orders,
    product_orders,
    reactant_symbols,
    product_indices_by_symbol,
):
    unassigned = order[pos:]
    if not unassigned:
        return 0.0

    remaining_by_symbol = {}
    for symbol, indices in product_indices_by_symbol.items():
        remaining_by_symbol[symbol] = [idx for idx in indices if not used_product[idx]]

    needed = Counter(reactant_symbols[idx] for idx in unassigned)
    for symbol, count in needed.items():
        if len(remaining_by_symbol[symbol]) < count:
            return None

    symbol_pair_minmax = {}
    min_total = 0.0
    assigned = order[:pos]
    minmax_assigned = {}

    for i in unassigned:
        sym_i = reactant_symbols[i]
        candidates_i = remaining_by_symbol[sym_i]
        if not candidates_i:
            return None
        for j in assigned:
            prod_j = mapping[j]
            cache_key = (prod_j, sym_i)
            if cache_key not in minmax_assigned:
                min_p, max_p = minmax_product_orders(
                    [prod_j], candidates_i, product_orders, allow_same=False
                )
                minmax_assigned[cache_key] = (min_p, max_p)
            min_p, max_p = minmax_assigned[cache_key]
            min_delta = min_p - reactant_orders[i][j]
            max_delta = max_p - reactant_orders[i][j]
            a0, a1, b0, b1 = add_break_bounds_from_delta(min_delta, max_delta)
            min_total += a0 + b0

    for idx_a, i in enumerate(unassigned):
        sym_i = reactant_symbols[i]
        for j in unassigned[idx_a + 1 :]:
            sym_j = reactant_symbols[j]
            key = (sym_i, sym_j)
            if key not in symbol_pair_minmax:
                allow_same = sym_i != sym_j
                min_p, max_p = minmax_product_orders(
                    remaining_by_symbol[sym_i],
                    remaining_by_symbol[sym_j],
                    product_orders,
                    allow_same=allow_same,
                )
                if min_p is None:
                    return None
                symbol_pair_minmax[key] = (min_p, max_p)
                symbol_pair_minmax[(sym_j, sym_i)] = (min_p, max_p)
            min_p, max_p = symbol_pair_minmax[key]
            min_delta = min_p - reactant_orders[i][j]
            max_delta = max_p - reactant_orders[i][j]
            a0, a1, b0, b1 = add_break_bounds_from_delta(min_delta, max_delta)
            min_total += a0 + b0

    return min_total


def build_mapped_reaction_smiles(reactant_mol, product_mol, atom_mapping):
    reactant = Chem.Mol(reactant_mol)
    product = Chem.Mol(product_mol)

    for atom in reactant.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    for start_idx, target_idx in atom_mapping.items():
        product.GetAtomWithIdx(target_idx).SetAtomMapNum(start_idx + 1)

    # Sanitize to get aromatic form output
    # Chem.SanitizeMol(reactant)
    # Chem.SanitizeMol(product)

    reactant_smiles = Chem.MolToSmiles(reactant, canonical=True, allHsExplicit=True)
    product_smiles = Chem.MolToSmiles(product, canonical=True, allHsExplicit=True)
    return f"{reactant_smiles}>>{product_smiles}"


def build_draw_path(base_path, index):
    base, ext = os.path.splitext(base_path)
    if not ext:
        ext = ".svg"
    if index == 0:
        return f"{base}{ext}"
    return f"{base}_{index + 1}{ext}"


def ensure_parent_dir(path):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)


def save_reaction_svg(reactant_mol, product_mol, atom_mapping, output_path):
    reactant = Chem.Mol(reactant_mol)
    product = Chem.Mol(product_mol)

    for atom in reactant.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    for start_idx, target_idx in atom_mapping.items():
        product.GetAtomWithIdx(target_idx).SetAtomMapNum(start_idx + 1)

    rxn = AllChem.ChemicalReaction()
    rxn.AddReactantTemplate(reactant)
    rxn.AddProductTemplate(product)
    svg = Draw.ReactionToImage(rxn, useSVG=True)
    ensure_parent_dir(output_path)
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(svg)


def _format_debug_done(nodes, kept, drawn, best_total, best_pairs, elapsed):
    best_total_str = "none" if best_total is None else format_count(best_total)
    best_pairs_str = format_best_pairs(best_pairs)
    return (
        f"[debug] done nodes={nodes} kept={kept} drawn={drawn} "
        f"best_total={best_total_str} best_ops={best_pairs_str} "
        f"elapsed={elapsed:.1f}s"
    )


def _apply_assignment(
    mapping, used_product, react_idx, prod_idx, reactant_orders, product_orders
):
    delta_add = 0.0
    delta_break = 0.0
    for assigned_idx, assigned_prod in mapping.items():
        delta = product_orders[prod_idx][assigned_prod] - reactant_orders[react_idx][
            assigned_idx
        ]
        if delta > 0:
            delta_add += delta
        elif delta < 0:
            delta_break += -delta

    mapping[react_idx] = prod_idx
    used_product[prod_idx] = True
    return delta_add, delta_break


def _run_search(
    reactant_mol,
    product_mol,
    *,
    order,
    mapping,
    used_product,
    pos0,
    add0,
    break0,
    reactant_orders,
    product_orders,
    reactant_symbols,
    product_indices_by_symbol,
    reactant_ranks,
    product_ranks,
    best_total0=None,
    best_pairs0=None,
    debug=False,
    debug_interval=100000,
    draw_path=None,
    draw_limit=0,
):
    results_by_key = {}

    nodes = 0
    kept = 0
    drawn = 0
    best_total = best_total0
    best_pairs = set(best_pairs0) if best_pairs0 else set()
    start_time = time.time()

    def backtrack(pos, add_count, break_count):
        nonlocal nodes, kept, drawn, best_total, best_pairs
        nodes += 1
        if debug and nodes % debug_interval == 0:
            elapsed = time.time() - start_time
            print(
                f"[debug] nodes={nodes} kept={kept} elapsed={elapsed:.1f}s",
                file=sys.stderr,
            )

        current_total = add_count + break_count
        if best_total is not None and current_total > best_total + EPS:
            return

        min_remaining = compute_remaining_min_total(
            order,
            pos,
            mapping,
            used_product,
            reactant_orders,
            product_orders,
            reactant_symbols,
            product_indices_by_symbol,
        )
        if min_remaining is None:
            return
        if best_total is not None and current_total + min_remaining > best_total + EPS:
            return

        if pos == len(order):
            total = current_total
            if best_total is None or total < best_total - EPS:
                best_total = total
                best_pairs = {(break_count, add_count)}
                results_by_key.clear()
            elif abs(total - best_total) <= EPS:
                best_pairs.add((break_count, add_count))
            else:
                return
            
            key = mapping_key(mapping, reactant_ranks, product_ranks)
            if key in results_by_key:
                return

            mapped = build_mapped_reaction_smiles(reactant_mol, product_mol, mapping)
            results_by_key[key] = mapped
            # if mapped in results_by_key:
                # return
            # results_by_key[mapped] = mapped
            kept += 1
            if draw_path and (draw_limit <= 0 or drawn < draw_limit):
                output_path = build_draw_path(draw_path, drawn)
                save_reaction_svg(reactant_mol, product_mol, mapping, output_path)
                drawn += 1
            return

        react_idx = order[pos]
        symbol = reactant_symbols[react_idx]
        candidates = product_indices_by_symbol[symbol]
        for prod_idx in candidates:
            if used_product[prod_idx]:
                continue
            delta_add, delta_break = _apply_assignment(
                mapping,
                used_product,
                react_idx,
                prod_idx,
                reactant_orders,
                product_orders,
            )
            backtrack(pos + 1, add_count + delta_add, break_count + delta_break)
            used_product[prod_idx] = False
            mapping.pop(react_idx, None)

    backtrack(pos0, add0, break0)
    elapsed = time.time() - start_time
    return {
        "results_by_key": results_by_key,
        "nodes": nodes,
        "kept": kept,
        "drawn": drawn,
        "best_total": best_total,
        "best_pairs": best_pairs,
        "elapsed": elapsed,
    }


def _find_first_branch(order, reactant_symbols, product_indices_by_symbol):
    for pos, react_idx in enumerate(order):
        symbol = reactant_symbols[react_idx]
        candidates = product_indices_by_symbol[symbol]
        if len(candidates) > 1:
            return pos, react_idx, candidates
    return None


def _worker_branch(args):
    rxn_smiles, chosen_prod_idx = args
    reactant_smiles, product_smiles = parse_reaction_smiles(rxn_smiles)
    reactant_mol = smiles_to_mol(reactant_smiles)
    product_mol = smiles_to_mol(product_smiles)

    if reactant_mol.GetNumAtoms() != product_mol.GetNumAtoms():
        raise ValueError("Reactant and product must have the same number of atoms.")
    if atom_symbol_counts(reactant_mol) != atom_symbol_counts(product_mol):
        raise ValueError("Reactant and product must have the same chemical formula.")

    reactant_graph = mol_to_graph(reactant_mol)
    reactant_symbols = [atom.GetSymbol() for atom in reactant_mol.GetAtoms()]
    product_symbols = [atom.GetSymbol() for atom in product_mol.GetAtoms()]
    product_indices_by_symbol = {}
    for idx, sym in enumerate(product_symbols):
        product_indices_by_symbol.setdefault(sym, []).append(idx)

    reactant_orders = bond_order_matrix(reactant_mol)
    product_orders = bond_order_matrix(product_mol)
    reactant_ranks = symmetry_ranks(reactant_mol)
    product_ranks = symmetry_ranks(product_mol)

    best_total0, best_pairs0 = _seed_best_from_input_mapping(
        reactant_mol,
        product_mol,
        reactant_orders=reactant_orders,
        product_orders=product_orders,
        reactant_symbols=reactant_symbols,
        product_indices_by_symbol=product_indices_by_symbol,
    )

    order = sorted(
        range(reactant_mol.GetNumAtoms()),
        key=lambda i: (
            len(product_indices_by_symbol[reactant_symbols[i]]),
            -reactant_graph.degree[i],
        ),
    )

    branch = _find_first_branch(order, reactant_symbols, product_indices_by_symbol)
    if branch is None:
        return {
            "results_by_key": {},
            "nodes": 0,
            "kept": 0,
            "drawn": 0,
            "best_total": None,
            "best_pairs": set(),
            "elapsed": 0.0,
        }

    split_pos, split_react_idx, candidates = branch
    candidates_set = set(candidates)
    if chosen_prod_idx not in candidates_set:
        raise ValueError("Invalid branch candidate for split position.")

    used_product = [False] * product_mol.GetNumAtoms()
    mapping = {}

    add_count = 0.0
    break_count = 0.0

    for pos in range(split_pos):
        react_idx = order[pos]
        symbol = reactant_symbols[react_idx]
        forced_candidates = product_indices_by_symbol[symbol]
        if len(forced_candidates) != 1:
            raise ValueError("Unexpected non-forced prefix while splitting.")
        prod_idx = forced_candidates[0]
        if used_product[prod_idx]:
            return {
                "results_by_key": {},
                "nodes": 0,
                "kept": 0,
                "drawn": 0,
                "best_total": None,
                "best_pairs": set(),
                "elapsed": 0.0,
            }
        da, db = _apply_assignment(
            mapping, used_product, react_idx, prod_idx, reactant_orders, product_orders
        )
        add_count += da
        break_count += db

    if used_product[chosen_prod_idx]:
        return {
            "results_by_key": {},
            "nodes": 0,
            "kept": 0,
            "drawn": 0,
            "best_total": None,
            "best_pairs": set(),
            "elapsed": 0.0,
        }
    da, db = _apply_assignment(
        mapping,
        used_product,
        split_react_idx,
        chosen_prod_idx,
        reactant_orders,
        product_orders,
    )
    add_count += da
    break_count += db

    return _run_search(
        reactant_mol,
        product_mol,
        order=order,
        mapping=mapping,
        used_product=used_product,
        pos0=split_pos + 1,
        add0=add_count,
        break0=break_count,
        reactant_orders=reactant_orders,
        product_orders=product_orders,
        reactant_symbols=reactant_symbols,
        product_indices_by_symbol=product_indices_by_symbol,
        reactant_ranks=reactant_ranks,
        product_ranks=product_ranks,
        best_total0=best_total0,
        best_pairs0=best_pairs0,
        debug=False,
        debug_interval=100000,
        draw_path=None,
        draw_limit=0,
    )


def generate_atom_mapped_reactions(
    rxn_smiles,
    debug=False,
    debug_interval=100000,
    draw_path=None,
    draw_limit=0,
    jobs=1,
):
    reactant_smiles, product_smiles = parse_reaction_smiles(rxn_smiles)
    reactant_mol = smiles_to_mol(reactant_smiles)
    product_mol = smiles_to_mol(product_smiles)

    # Check atom number and formula
    if reactant_mol.GetNumAtoms() != product_mol.GetNumAtoms():
        raise ValueError("Reactant and product must have the same number of atoms.")
    if atom_symbol_counts(reactant_mol) != atom_symbol_counts(product_mol):
        raise ValueError("Reactant and product must have the same chemical formula.")

    reactant_graph = mol_to_graph(reactant_mol)
    product_graph = mol_to_graph(product_mol)
    reactant_groups = group_atoms_by_symbol(reactant_graph)
    product_groups = group_atoms_by_symbol(product_graph)
    if reactant_groups.keys() != product_groups.keys():
        return []

    reactant_symbols = [atom.GetSymbol() for atom in reactant_mol.GetAtoms()]
    product_symbols = [atom.GetSymbol() for atom in product_mol.GetAtoms()]
    product_indices_by_symbol = {}
    for idx, sym in enumerate(product_symbols):
        product_indices_by_symbol.setdefault(sym, []).append(idx)

    reactant_orders = bond_order_matrix(reactant_mol)
    product_orders = bond_order_matrix(product_mol)
    reactant_ranks = symmetry_ranks(reactant_mol)
    product_ranks = symmetry_ranks(product_mol)

    best_total0, best_pairs0 = _seed_best_from_input_mapping(
        reactant_mol,
        product_mol,
        reactant_orders=reactant_orders,
        product_orders=product_orders,
        reactant_symbols=reactant_symbols,
        product_indices_by_symbol=product_indices_by_symbol,
    )
    if debug and best_total0 is not None:
        best_total_str = format_count(best_total0)
        best_pairs_str = format_best_pairs(best_pairs0)
        print(
            f"[debug] seed best_total={best_total_str} best_ops={best_pairs_str}",
            file=sys.stderr,
        )
    elif debug:
        reason = _describe_seed_skip_reason(reactant_mol, product_mol)
        print(f"[debug] seed skipped: {reason}", file=sys.stderr)

    if jobs == 0:
        jobs = os.cpu_count() or 1
    if jobs < 1:
        raise ValueError("--jobs must be >= 1 (or 0 for all CPUs).")
    if jobs > 1 and draw_path:
        raise ValueError("Drawing is not supported with --jobs > 1.")

    if debug:
        counts = atom_symbol_counts(reactant_mol)
        group_sizes = {symbol: len(nodes) for symbol, nodes in reactant_groups.items()}
        est = estimate_mapping_count(reactant_groups)
        print(f"[debug] atoms: {reactant_mol.GetNumAtoms()}", file=sys.stderr)
        print(f"[debug] symbol_counts: {dict(counts)}", file=sys.stderr)
        print(f"[debug] group_sizes: {group_sizes}", file=sys.stderr)
        print(f"[debug] estimated_mappings: {est}", file=sys.stderr)

    order = sorted(
        range(reactant_mol.GetNumAtoms()),
        key=lambda i: (
            len(product_indices_by_symbol[reactant_symbols[i]]),
            -reactant_graph.degree[i],
        ),
    )

    if jobs == 1:
        used_product = [False] * product_mol.GetNumAtoms()
        mapping = {}
        out = _run_search(
            reactant_mol,
            product_mol,
            order=order,
            mapping=mapping,
            used_product=used_product,
            pos0=0,
            add0=0.0,
            break0=0.0,
            reactant_orders=reactant_orders,
            product_orders=product_orders,
            reactant_symbols=reactant_symbols,
            product_indices_by_symbol=product_indices_by_symbol,
            reactant_ranks=reactant_ranks,
            product_ranks=product_ranks,
            best_total0=best_total0,
            best_pairs0=best_pairs0,
            debug=debug,
            debug_interval=debug_interval,
            draw_path=draw_path,
            draw_limit=draw_limit,
        )
        if debug:
            print(
                _format_debug_done(
                    out["nodes"],
                    out["kept"],
                    out["drawn"],
                    out["best_total"],
                    out["best_pairs"],
                    out["elapsed"],
                ),
                file=sys.stderr,
            )
        return sorted(out["results_by_key"].values())

    branch = _find_first_branch(order, reactant_symbols, product_indices_by_symbol)
    if branch is None:
        used_product = [False] * product_mol.GetNumAtoms()
        mapping = {}
        out = _run_search(
            reactant_mol,
            product_mol,
            order=order,
            mapping=mapping,
            used_product=used_product,
            pos0=0,
            add0=0.0,
            break0=0.0,
            reactant_orders=reactant_orders,
            product_orders=product_orders,
            reactant_symbols=reactant_symbols,
            product_indices_by_symbol=product_indices_by_symbol,
            reactant_ranks=reactant_ranks,
            product_ranks=product_ranks,
            best_total0=best_total0,
            best_pairs0=best_pairs0,
            debug=debug,
            debug_interval=debug_interval,
            draw_path=draw_path,
            draw_limit=draw_limit,
        )
        if debug:
            print(
                _format_debug_done(
                    out["nodes"],
                    out["kept"],
                    out["drawn"],
                    out["best_total"],
                    out["best_pairs"],
                    out["elapsed"],
                ),
                file=sys.stderr,
            )
        return sorted(out["results_by_key"].values())

    split_pos, split_react_idx, candidates = branch
    if debug:
        symbol = reactant_symbols[split_react_idx]
        max_workers = min(jobs, len(candidates))
        print(
            f"[debug] parallel split pos={split_pos} react_idx={split_react_idx} "
            f"symbol={symbol} branches={len(candidates)} jobs={jobs} workers={max_workers}",
            file=sys.stderr,
        )

    tasks = [(rxn_smiles, prod_idx) for prod_idx in candidates]
    max_workers = min(jobs, len(tasks))
    merged_by_key = {}
    merged_nodes = 0
    best_total = None
    best_pairs = set()
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(_worker_branch, task) for task in tasks]
        for fut in as_completed(futures):
            out = fut.result()
            merged_nodes += out["nodes"]
            local_best = out["best_total"]
            if local_best is None:
                continue
            if best_total is None or local_best < best_total - EPS:
                best_total = local_best
                best_pairs = set(out["best_pairs"])
                merged_by_key.clear()
                merged_by_key.update(out["results_by_key"])
            elif abs(local_best - best_total) <= EPS:
                best_pairs |= set(out["best_pairs"])
                for key, mapped in out["results_by_key"].items():
                    merged_by_key.setdefault(key, mapped)

    elapsed = time.time() - start_time
    if debug:
        print(
            _format_debug_done(
                merged_nodes,
                len(merged_by_key),
                0,
                best_total,
                best_pairs,
                elapsed,
            ),
            file=sys.stderr,
        )
    return sorted(merged_by_key.values())


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate unique atom-mapped reaction SMILES with minimal add+break count."
        )
    )
    parser.add_argument("smiles", help="Reaction SMILES, e.g. O=CCCOO>>CC=O.O=CO")
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug information and progress to stderr.",
    )
    parser.add_argument(
        "--debug-interval",
        type=int,
        default=100000,
        help="Progress interval when --debug is set, default 100000.",
    )
    parser.add_argument(
        "--draw",
        action="store_true",
        help="Save reaction SVGs for accepted mappings.",
    )
    parser.add_argument(
        "--draw-path",
        default="reaction.svg",
        help="Output SVG path or base name when --draw is set.",
    )
    parser.add_argument(
        "--draw-limit",
        type=int,
        default=0,
        help="Maximum number of SVGs to write when --draw is set; 0 for all.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of worker processes for a single reaction; 0 uses all CPUs.",
    )
    args = parser.parse_args()

    results = generate_atom_mapped_reactions(
        args.smiles,
        debug=args.debug,
        debug_interval=args.debug_interval,
        draw_path=args.draw_path if args.draw else None,
        draw_limit=args.draw_limit,
        jobs=args.jobs,
    )

    for rxn in results:
        print(rxn)


if __name__ == "__main__":
    main()
