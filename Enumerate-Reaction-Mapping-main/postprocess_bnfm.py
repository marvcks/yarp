"""
Post-process BNFM enumeration results.

Steps:
1. Flatten reactant -> products into reactant/product pairs.
2. Deduplicate reactions by unordered (reactant_hash, product_hash) so A->B and B->A are treated as duplicates.
3. Remove `geo` from all species records.
4. Generate plain and atom-mapped reaction SMILES (RDKit-based).
5. Save cleaned reactions to a pickle and print unique species/reaction counts.
"""

from __future__ import annotations

import argparse
import array
from concurrent.futures import ProcessPoolExecutor
import hashlib
import logging
import multiprocessing as mp
import os
import pickle
import tempfile
import time
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
from rdkit import Chem
import yarp as yp
from yarp.find_lewis import return_formals

try:
    from tqdm import tqdm  # type: ignore
except Exception:  # pragma: no cover
    tqdm = None

try:
    import resource  # type: ignore
except Exception:  # pragma: no cover
    resource = None  # type: ignore

logger = logging.getLogger(__name__)


def return_smi(elements, bond_mat, fc):
    """
    Generate atom-mapped and plain SMILES from elements/bond_mat/fc.
    """
    try:
        mol = Chem.RWMol()

        # Add atoms with mapping numbers
        for i, element in enumerate(elements):
            atom = Chem.Atom(element)
            atom.SetAtomMapNum(i + 1)
            atom.SetFormalCharge(int(fc[i]))
            atom.SetNumRadicalElectrons(int(bond_mat[i, i] % 2))
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
            mol.AddAtom(atom)

        # Add bonds
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
                # No other bond types handled
                # If other types (e.g., aromatic, quadruple) are needed, add here
                else:
                    continue
                mol.AddBond(i, j, bond_type)

        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol)

        mapped = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)

        mol_plain = Chem.RemoveHs(mol)
        for atom in mol_plain.GetAtoms():
            atom.SetAtomMapNum(0)
        plain = Chem.MolToSmiles(mol_plain, canonical=True, isomericSmiles=False)

        return mapped, plain
    except Exception as exc:
        print(f"Error generating SMILES: {exc}")
        return "", ""


def strip_geo(spec: Dict) -> Dict:
    out = dict(spec)
    out.pop("geo", None)
    return out


def reaction_key(h1: float, h2: float) -> Tuple[float, float]:
    return (h1, h2) if h1 <= h2 else (h2, h1)


def _maxrss_mb() -> Optional[float]:
    if resource is None:
        return None
    try:
        maxrss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # Linux: kilobytes; macOS: bytes. Use a simple heuristic.
        if maxrss > 100_000_000:
            return maxrss / (1024 * 1024)
        return maxrss / 1024.0
    except Exception:
        return None


def _union_len(a: Dict[float, object], b: Dict[float, object]) -> int:
    if not a:
        return len(b)
    if not b:
        return len(a)
    if len(a) <= len(b):
        small, large = a, b
    else:
        small, large = b, a
    inter = 0
    for k in small.keys():
        if k in large:
            inter += 1
    return len(a) + len(b) - inter


def _spec_shape_sig(spec: Dict) -> Tuple[int, Tuple[int, ...]]:
    """
    Very cheap signature to detect obvious hash collisions / record mismatches:
    (#atoms, bond_mats shape).
    """
    n_atoms = len(spec.get("elements", []))
    try:
        shape = tuple(np.asarray(spec.get("bond_mats")).shape)
    except Exception:
        shape = tuple()
    return n_atoms, shape


def _spec_fingerprint_u64(spec: Dict) -> int:
    """
    Order-sensitive fingerprint for a species record (elements + bond_mats raw bytes).
    Used for variant caching and debugging/logging.
    """
    h = hashlib.blake2b(digest_size=8)
    elems = spec.get("elements", [])
    h.update(str(tuple(elems)).encode("utf-8", errors="ignore"))
    try:
        b = np.asarray(spec.get("bond_mats"))
        h.update(str(b.shape).encode("ascii"))
        h.update(str(b.dtype).encode("ascii"))
        b_contig = np.ascontiguousarray(b)
        h.update(memoryview(b_contig))
    except Exception:
        pass
    return int.from_bytes(h.digest(), byteorder="little", signed=False)


def _bond_mat_from_record(bond_mats) -> np.ndarray:
    bmat_raw = np.asarray(bond_mats)
    if bmat_raw.ndim == 3:
        return bmat_raw[0]
    elif bmat_raw.ndim == 2:
        return bmat_raw
    elif bmat_raw.ndim == 1:
        side = int(round(len(bmat_raw) ** 0.5))
        if side * side != len(bmat_raw):
            raise ValueError(f"Unable to reshape bond_mats of length {len(bmat_raw)} into square matrix.")
        return bmat_raw.reshape((side, side))
    else:
        raise ValueError(f"Unexpected bond_mats shape: {bmat_raw.shape}")


def species_smiles(spec: Dict) -> Tuple[str, str]:
    elements = [element.upper() for element in spec["elements"]]
    bmat = _bond_mat_from_record(spec["bond_mats"])
    fc = return_formals(bmat, elements)
    return return_smi(elements, bmat, fc)


def _species_smiles_from_parts(elements: Sequence[str], bond_mats) -> Tuple[str, str]:
    elements_up = [element.upper() for element in elements]
    bmat = _bond_mat_from_record(bond_mats)
    fc = return_formals(bmat, elements_up)
    return return_smi(elements_up, bmat, fc)


def _iter_bnfm_records(file_obj) -> Iterator[Dict]:
    """
    Yield BNFM reactant records from either:
    - a single pickled list
    - or a stream of pickled lists (one list per iteration)
    """
    while True:
        try:
            chunk = pickle.load(file_obj)
        except EOFError:
            break
        if isinstance(chunk, list):
            for entry in chunk:
                yield entry
        else:
            yield chunk


def _load_seed_smiles_file(path: str) -> List[str]:
    smiles: List[str] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smiles.extend(line.split())
    return smiles


def _load_seed_smiles_by_formula(path: str, formula: str) -> List[str]:
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts and parts[0] == formula:
                return parts[1:]
    raise ValueError(f"Formula {formula!r} not found in {path!r}")


def _load_seed_smiles_by_line(path: str, line_no: int) -> Tuple[str, List[str]]:
    """
    Load seed SMILES from a smiles_by_formula.txt-style file by 1-based line number.
    Returns (formula, smiles_list).
    """
    if line_no <= 0:
        raise ValueError("line_no must be >= 1")
    with open(path) as f:
        for idx, line in enumerate(f, start=1):
            if idx != line_no:
                continue
            line = line.strip()
            parts = line.split()
            if not parts:
                raise ValueError(f"Line {line_no} is empty in {path!r}")
            return parts[0], parts[1:]
    raise ValueError(f"Line {line_no} not found in {path!r}")


def _hashes_from_smiles(smiles: Sequence[str]) -> set[float]:
    hashes: set[float] = set()
    for smi in smiles:
        smi = smi.strip()
        if not smi:
            continue
        try:
            mol = yp.yarpecule(smi)
            hashes.add(float(mol.hash))
        except Exception as exc:
            logger.warning("Seed SMILES parse failed: %s | %s", smi, exc)
    return hashes


def _species_smiles_task(task: Tuple[float, Sequence[str], object]) -> Tuple[float, str, str, str]:
    """
    Compute SMILES for a single species record.

    Returns (hash, mapped, plain, error_message). error_message is empty string when successful.
    """
    species_hash, elements, bond_mats = task
    try:
        mapped, plain = _species_smiles_from_parts(elements, bond_mats)
        return species_hash, mapped, plain, ""
    except Exception as exc:
        return species_hash, "", "", str(exc)


def _dedup_bucket_worker(
    bucket_id: int,
    in_queue: "mp.queues.Queue[Optional[bytes]]",
    out_queue: "mp.queues.Queue[Tuple[int, str, int]]",
    tmp_dir: str,
    with_variant_ids: bool = False,
) -> None:
    """
    Worker that deduplicates unordered reaction pairs for a single bucket.

    Receives batches as raw float64 bytes:
    - with_variant_ids=False: (r_hash, p_hash, r_hash, p_hash, ...)
    - with_variant_ids=True:  (r_hash, p_hash, r_var_id, p_var_id, ...)

    Writes:
    - with_variant_ids=False: float64 pairs (r_hash, p_hash) for each unique unordered key
    - with_variant_ids=True:  uint64 pairs (r_var_id, p_var_id) for each unique unordered key

    Returns (bucket_id, path, num_pairs).
    """
    seen: set[Tuple[float, float]] = set()

    out_path = os.path.join(tmp_dir, f"dedup_bucket_{bucket_id:04d}.bin")
    pairs_written = 0
    # Write incrementally to keep worker RSS bounded.
    with open(out_path, "wb") as f:
        if with_variant_ids:
            buf_ids = array.array("Q")
        else:
            buf = array.array("d")
        while True:
            payload = in_queue.get()
            if payload is None:
                break
            arr = np.frombuffer(payload, dtype=np.float64)
            if with_variant_ids:
                if arr.size % 4 != 0:
                    raise ValueError(f"Unexpected payload length (not quads): {arr.size}")
                # Iterate quads: arr is [r0, p0, rid0, pid0, r1, p1, rid1, pid1, ...]
                for i in range(0, arr.size, 4):
                    r_hash = arr[i]
                    p_hash = arr[i + 1]
                    key = (r_hash, p_hash) if r_hash <= p_hash else (p_hash, r_hash)
                    if key in seen:
                        continue
                    seen.add(key)
                    buf_ids.append(int(arr[i + 2]))
                    buf_ids.append(int(arr[i + 3]))

                # Flush occasionally so the buffer doesn't grow without bound.
                if buf_ids and len(buf_ids) >= 2_000_000:
                    f.write(buf_ids.tobytes())
                    pairs_written += len(buf_ids) // 2
                    buf_ids = array.array("Q")
            else:
                if arr.size % 2 != 0:
                    raise ValueError(f"Unexpected payload length (not pairs): {arr.size}")
                # Iterate pairs: arr is [r0, p0, r1, p1, ...]
                for i in range(0, arr.size, 2):
                    r_hash = arr[i]
                    p_hash = arr[i + 1]
                    key = (r_hash, p_hash) if r_hash <= p_hash else (p_hash, r_hash)
                    if key in seen:
                        continue
                    seen.add(key)
                    buf.append(float(r_hash))
                    buf.append(float(p_hash))

                # Flush occasionally so the buffer doesn't grow without bound.
                if buf and len(buf) >= 2_000_000:
                    f.write(buf.tobytes())
                    pairs_written += len(buf) // 2
                    buf = array.array("d")

        if with_variant_ids:
            if buf_ids:
                f.write(buf_ids.tobytes())
                pairs_written += len(buf_ids) // 2
        else:
            if buf:
                f.write(buf.tobytes())
                pairs_written += len(buf) // 2

    out_queue.put((bucket_id, out_path, pairs_written))


def process(
    input_path: str,
    output_path: str,
    use_tqdm: bool = True,
    max_workers: int = 0,
    chunksize: int = 50,
    log_every_pairs: int = 1_000_000,
    dedup_workers: int = 1,
    dedup_batch_pairs: int = 50_000,
    max_depth: int = -1,
    depth0_seed_hashes: Optional[set[float]] = None,
    strict_role_records: bool = False,
    verify_output_len: bool = False,
    smiles_from_entry: bool = False,
) -> None:
    disable_tqdm = (not use_tqdm) or (tqdm is None)
    run_start = time.time()
    smiles_elapsed = 0.0
    build_elapsed = 0.0

    logger.info("Reading input pkl: %s", input_path)
    depth_filter = max_depth >= 0
    seed_candidates: set[float] = set(depth0_seed_hashes or set())
    if depth_filter and not seed_candidates:
        raise ValueError("--max-depth requires non-empty depth0_seed_hashes (seed SMILES).")
    if depth_filter:
        logger.info(
            "Depth filter enabled | max_depth=%d | seed_candidates=%d",
            max_depth,
            len(seed_candidates),
        )
    if smiles_from_entry:
        if strict_role_records:
            logger.warning("--strict-role-records is ignored when --smiles-from-entry is enabled.")
        logger.info(
            "SMILES mode: from entry records (variant-aware). "
            "This is slower to build variant cache, but preserves atom-map correspondence per r->p entry."
        )

    # Keep memory low: store one record per species hash, split by role.
    # Atom-mapped SMILES relies on consistent atom ordering; reactant/product role can differ in stored pkl.
    reactant_by_hash: Dict[float, Dict] = {}
    product_by_hash: Dict[float, Dict] = {}
    # Debug/guardrail: detect obvious per-hash record mismatches (e.g., float-hash collision).
    reactant_shape_sig_by_hash: Dict[float, Tuple[int, Tuple[int, ...]]] = {}
    product_shape_sig_by_hash: Dict[float, Tuple[int, Tuple[int, ...]]] = {}
    # Variant cache (for --smiles-from-entry): do NOT assume hash->record is unique or stable across depths.
    variant_id_by_fp: Dict[int, int] = {}
    variant_specs: List[Optional[Dict]] = [None]  # 1-based ids; slot 0 unused
    dedup_reactions: List[Tuple[float, float]] = []
    seen_reaction_keys: set[Tuple[float, float]] = set()
    reactant_records = 0
    reaction_pairs = 0

    depth = 0
    depth_start_t = time.time()
    depth0_reactants: set[float] = set()
    depth0_remaining: Optional[set[float]] = set(seed_candidates) if depth_filter else None
    # Treat all seed candidates as "already known species" so that depth-0 reactants discovered later
    # won't be incorrectly promoted into depth-1 when they appear as products earlier in the file.
    seen_species: set[float] = set(seed_candidates) if depth_filter else set()
    next_frontier: set[float] = set()
    frontier_remaining: Optional[set[float]] = None
    frontier_size = 0
    if depth_filter:
        logger.info("Depth 0 start | seed_candidates=%d", len(seed_candidates))

    if dedup_workers > 1:
        logger.info(
            "Dedup mode: bucketed multiprocessing | workers=%d | batch_pairs=%d",
            dedup_workers,
            dedup_batch_pairs,
        )
        ctx = mp.get_context("spawn")
        with tempfile.TemporaryDirectory(prefix="bnfm_dedup_") as tmp_dir:
            in_queues: List[mp.Queue[Optional[bytes]]] = [
                ctx.Queue(maxsize=4) for _ in range(dedup_workers)
            ]
            out_queue: mp.Queue[Tuple[int, str, int]] = ctx.Queue()
            procs: List[mp.Process] = []
            for bucket_id in range(dedup_workers):
                p = ctx.Process(
                    target=_dedup_bucket_worker,
                    args=(bucket_id, in_queues[bucket_id], out_queue, tmp_dir, smiles_from_entry),
                    daemon=True,
                )
                p.start()
                procs.append(p)

            # Per-bucket outgoing buffers: raw float64 pairs, to minimize IPC overhead.
            bucket_buffers: List[array.array] = [array.array("d") for _ in range(dedup_workers)]

            with open(input_path, "rb") as f:
                read_start = time.time()
                last_log_t = read_start
                last_log_pairs = 0
                next_log_pairs = log_every_pairs if log_every_pairs > 0 else None

                reactant_by_hash_local = reactant_by_hash
                product_by_hash_local = product_by_hash
                reactant_shape_sig_local = reactant_shape_sig_by_hash
                product_shape_sig_local = product_shape_sig_by_hash
                variant_id_by_fp_local = variant_id_by_fp
                variant_specs_local = variant_specs

                records2: Iterable[Dict] = _iter_bnfm_records(f)
                if tqdm is not None:
                    records2 = tqdm(records2, desc="Reading/Dedup", unit="reactant", disable=disable_tqdm)
                stop_reading = False
                for entry in records2:
                    reactant_records += 1
                    reactant_raw = entry["reactant"]
                    r_hash = reactant_raw["hash"]

                    # Depth filtering: treat depth 0 as "reactants whose hash is in seed_candidates".
                    if depth_filter:
                        while True:
                            if depth == 0:
                                if depth0_remaining is None:
                                    raise RuntimeError("Internal error: depth0_remaining is None while depth==0.")
                                if r_hash not in seed_candidates:
                                    # If we haven't seen any seed reactants yet, skip until we hit the first seed.
                                    if len(depth0_remaining) == len(seed_candidates):
                                        continue
                                    seen_cnt = len(seed_candidates) - len(depth0_remaining)
                                    raise ValueError(
                                        "Depth filter mismatch: encountered a non-seed reactant before consuming all seed reactants. "
                                        f"seen_seed={seen_cnt}/{len(seed_candidates)} | remaining_seed={len(depth0_remaining)} | reactant_hash={r_hash}. "
                                        "This usually means the seed SMILES (formula/line) does not match this pkl, or the pkl contains mixed/resumed iterations."
                                    )
                                depth0_reactants.add(r_hash)
                                depth0_remaining.discard(r_hash)
                                break

                            # depth >= 1
                            if frontier_remaining is None:
                                raise RuntimeError("Internal error: frontier_remaining is None for depth>=1.")
                            if r_hash not in frontier_remaining:
                                raise ValueError(
                                    f"Depth filter mismatch at depth={depth}: reactant hash {r_hash} not in frontier."
                                )
                            frontier_remaining.remove(r_hash)
                            break
                        if stop_reading:
                            break

                    reactant_raw.pop("geo", None)
                    r_var_id: Optional[int] = None
                    if smiles_from_entry:
                        r_fp = _spec_fingerprint_u64(reactant_raw)
                        r_var_id = variant_id_by_fp_local.get(r_fp)
                        if r_var_id is None:
                            r_var_id = len(variant_specs_local)
                            variant_id_by_fp_local[r_fp] = r_var_id
                            variant_specs_local.append(reactant_raw)

                    if r_hash not in reactant_by_hash_local:
                        reactant_by_hash_local[r_hash] = reactant_raw
                        reactant_shape_sig_local[r_hash] = _spec_shape_sig(reactant_raw)
                    else:
                        prev_sig = reactant_shape_sig_local.get(r_hash)
                        cur_sig = _spec_shape_sig(reactant_raw)
                        if prev_sig is not None and prev_sig != cur_sig:
                            logger.warning(
                                "Reactant record mismatch for hash=%s | prev=%s cur=%s | prev_fp=%016x cur_fp=%016x",
                                r_hash,
                                prev_sig,
                                cur_sig,
                                _spec_fingerprint_u64(reactant_by_hash_local[r_hash]),
                                _spec_fingerprint_u64(reactant_raw),
                            )

                    products = entry.get("products")
                    if products:
                        for product_raw in products:
                            reaction_pairs += 1
                            p_hash = product_raw["hash"]
                            product_raw.pop("geo", None)
                            p_var_id: Optional[int] = None
                            if smiles_from_entry:
                                p_fp = _spec_fingerprint_u64(product_raw)
                                p_var_id = variant_id_by_fp_local.get(p_fp)
                                if p_var_id is None:
                                    p_var_id = len(variant_specs_local)
                                    variant_id_by_fp_local[p_fp] = p_var_id
                                    variant_specs_local.append(product_raw)

                            if p_hash not in product_by_hash_local:
                                product_by_hash_local[p_hash] = product_raw
                                product_shape_sig_local[p_hash] = _spec_shape_sig(product_raw)
                            else:
                                prev_sig = product_shape_sig_local.get(p_hash)
                                cur_sig = _spec_shape_sig(product_raw)
                                if prev_sig is not None and prev_sig != cur_sig:
                                    logger.warning(
                                        "Product record mismatch for hash=%s | prev=%s cur=%s | prev_fp=%016x cur_fp=%016x",
                                        p_hash,
                                        prev_sig,
                                        cur_sig,
                                        _spec_fingerprint_u64(product_by_hash_local[p_hash]),
                                        _spec_fingerprint_u64(product_raw),
                                    )

                            if depth_filter and p_hash not in seen_species:
                                seen_species.add(p_hash)
                                next_frontier.add(p_hash)

                            key = (r_hash, p_hash) if r_hash <= p_hash else (p_hash, r_hash)
                            bucket = hash(key) % dedup_workers
                            buf = bucket_buffers[bucket]
                            buf.append(r_hash)
                            buf.append(p_hash)
                            if smiles_from_entry:
                                if r_var_id is None or p_var_id is None:
                                    raise RuntimeError("Internal error: missing variant ids while smiles_from_entry.")
                                buf.append(float(r_var_id))
                                buf.append(float(p_var_id))

                            stride = 4 if smiles_from_entry else 2
                            if len(buf) >= stride * dedup_batch_pairs:
                                in_queues[bucket].put(buf.tobytes())
                                bucket_buffers[bucket] = array.array("d")

                            if next_log_pairs is not None and reaction_pairs >= next_log_pairs:
                                now = time.time()
                                dt = now - last_log_t
                                dp = reaction_pairs - last_log_pairs
                                pairs_per_s = (dp / dt) if dt > 0 else 0.0
                                total_dt = now - read_start
                                total_pairs_per_s = (reaction_pairs / total_dt) if total_dt > 0 else 0.0
                                rss_mb = _maxrss_mb()
                                rss_str = f"{rss_mb:.1f} MB" if rss_mb is not None else "n/a"
                                depth_str = f"{depth}" if depth_filter else "NA"
                                rem_str = "NA"
                                if depth_filter and depth >= 1 and frontier_remaining is not None:
                                    rem_str = str(len(frontier_remaining))
                                unique_species = _union_len(reactant_by_hash_local, product_by_hash_local)
                                variants_str = ""
                                if smiles_from_entry:
                                    variants_str = f" | variants={len(variant_specs_local) - 1}"
                                logger.info(
                                    "Reading/Dedup (mp) | depth=%s | pairs=%d | rate=%.0f pairs/s (avg %.0f) | frontier_rem=%s | unique_species=%d%s | rss=%s",
                                    depth_str,
                                    reaction_pairs,
                                    pairs_per_s,
                                    total_pairs_per_s,
                                    rem_str,
                                    unique_species,
                                    variants_str,
                                    rss_str,
                                )
                                last_log_t = now
                                last_log_pairs = reaction_pairs
                                next_log_pairs += log_every_pairs

                    if depth_filter and depth == 0 and depth0_remaining is not None and not depth0_remaining:
                        depth_elapsed = time.time() - depth_start_t
                        logger.info(
                            "Depth 0 done | reactants=%d | new_species=%d | elapsed %.2fs",
                            len(depth0_reactants),
                            len(next_frontier),
                            depth_elapsed,
                        )
                        if max_depth <= 0:
                            stop_reading = True
                            break
                        depth = 1
                        frontier_remaining = next_frontier
                        next_frontier = set()
                        frontier_size = len(frontier_remaining)
                        depth_start_t = time.time()
                        logger.info("Depth %d start | frontier=%d", depth, frontier_size)
                        if frontier_size == 0:
                            stop_reading = True
                            break

                    if depth_filter and depth >= 1 and frontier_remaining is not None and not frontier_remaining:
                        depth_elapsed = time.time() - depth_start_t
                        logger.info(
                            "Depth %d done | reactants=%d | new_species=%d | elapsed %.2fs",
                            depth,
                            frontier_size,
                            len(next_frontier),
                            depth_elapsed,
                        )
                        if depth >= max_depth:
                            stop_reading = True
                            break
                        depth += 1
                        frontier_remaining = next_frontier
                        next_frontier = set()
                        frontier_size = len(frontier_remaining)
                        depth_start_t = time.time()
                        logger.info("Depth %d start | frontier=%d", depth, frontier_size)
                        if frontier_size == 0:
                            stop_reading = True
                            break

                if stop_reading:
                    logger.info("Stopped reading early due to max_depth=%d", max_depth)

            # Flush remaining buffers.
            for bucket_id, buf in enumerate(bucket_buffers):
                if buf:
                    in_queues[bucket_id].put(buf.tobytes())
                    bucket_buffers[bucket_id] = array.array("d")

            # Tell workers to stop.
            for q in in_queues:
                q.put(None)

            # Collect worker outputs.
            bucket_outputs: List[Tuple[int, str, int]] = []
            for _ in range(dedup_workers):
                bucket_outputs.append(out_queue.get())
            bucket_outputs.sort(key=lambda t: t[0])

            # Join processes.
            for p in procs:
                p.join()

            # Load dedup pairs back into memory as a compact array.
            # Order is deterministic by bucket id, but not the original global-first-seen order.
            total_unique = sum(n_pairs for _, _, n_pairs in bucket_outputs)
            logger.info("Dedup done (mp) | unique_reactions=%d | loading buckets into memory", total_unique)
            dedup_dtype = np.uint64 if smiles_from_entry else np.float64
            dedup_pairs_arr = np.empty((total_unique, 2), dtype=dedup_dtype)
            offset = 0
            for bucket_id, path, n_pairs in bucket_outputs:
                if n_pairs == 0:
                    continue
                expected = 2 * n_pairs
                arr = np.fromfile(path, dtype=dedup_dtype, count=expected)
                if arr.size != expected:
                    kind = "uint64" if smiles_from_entry else "float64"
                    raise ValueError(f"Bucket {bucket_id}: expected {expected} {kind} values, got {arr.size} in {path}")
                dedup_pairs_arr[offset : offset + n_pairs, :] = arr.reshape((n_pairs, 2))
                offset += n_pairs
            if offset != total_unique:
                raise RuntimeError(f"Internal error: loaded {offset} pairs but expected {total_unique}")
            # Keep using the same variable name; numpy arrays are iterable as (r,p) rows.
            dedup_reactions = dedup_pairs_arr  # type: ignore[assignment]

    else:
        with open(input_path, "rb") as f:
            read_start = time.time()
            last_log_t = read_start
            last_log_pairs = 0
            next_log_pairs = log_every_pairs if log_every_pairs > 0 else None

            reactant_by_hash_local = reactant_by_hash
            product_by_hash_local = product_by_hash
            reactant_shape_sig_local = reactant_shape_sig_by_hash
            product_shape_sig_local = product_shape_sig_by_hash
            variant_id_by_fp_local = variant_id_by_fp
            variant_specs_local = variant_specs
            seen_reaction_keys_local = seen_reaction_keys
            seen_add = seen_reaction_keys_local.add
            dedup_append = dedup_reactions.append

            records: Iterable[Dict] = _iter_bnfm_records(f)
            if tqdm is not None:
                records = tqdm(records, desc="Reading/Dedup", unit="reactant", disable=disable_tqdm)
            stop_reading = False
            for entry in records:
                reactant_records += 1
                reactant_raw = entry["reactant"]
                r_hash = reactant_raw["hash"]

                if depth_filter:
                    while True:
                        if depth == 0:
                            if depth0_remaining is None:
                                raise RuntimeError("Internal error: depth0_remaining is None while depth==0.")
                            if r_hash not in seed_candidates:
                                if len(depth0_remaining) == len(seed_candidates):
                                    continue
                                seen_cnt = len(seed_candidates) - len(depth0_remaining)
                                raise ValueError(
                                    "Depth filter mismatch: encountered a non-seed reactant before consuming all seed reactants. "
                                    f"seen_seed={seen_cnt}/{len(seed_candidates)} | remaining_seed={len(depth0_remaining)} | reactant_hash={r_hash}. "
                                    "This usually means the seed SMILES (formula/line) does not match this pkl, or the pkl contains mixed/resumed iterations."
                                )
                            depth0_reactants.add(r_hash)
                            depth0_remaining.discard(r_hash)
                            break

                        if frontier_remaining is None:
                            raise RuntimeError("Internal error: frontier_remaining is None for depth>=1.")
                        if r_hash not in frontier_remaining:
                            raise ValueError(
                                f"Depth filter mismatch at depth={depth}: reactant hash {r_hash} not in frontier."
                            )
                        frontier_remaining.remove(r_hash)
                        break
                    if stop_reading:
                        break

                reactant_raw.pop("geo", None)
                r_var_id: Optional[int] = None
                if smiles_from_entry:
                    r_fp = _spec_fingerprint_u64(reactant_raw)
                    r_var_id = variant_id_by_fp_local.get(r_fp)
                    if r_var_id is None:
                        r_var_id = len(variant_specs_local)
                        variant_id_by_fp_local[r_fp] = r_var_id
                        variant_specs_local.append(reactant_raw)

                if r_hash not in reactant_by_hash_local:
                    reactant_by_hash_local[r_hash] = reactant_raw
                    reactant_shape_sig_local[r_hash] = _spec_shape_sig(reactant_raw)
                else:
                    prev_sig = reactant_shape_sig_local.get(r_hash)
                    cur_sig = _spec_shape_sig(reactant_raw)
                    if prev_sig is not None and prev_sig != cur_sig:
                        logger.warning(
                            "Reactant record mismatch for hash=%s | prev=%s cur=%s | prev_fp=%016x cur_fp=%016x",
                            r_hash,
                            prev_sig,
                            cur_sig,
                            _spec_fingerprint_u64(reactant_by_hash_local[r_hash]),
                            _spec_fingerprint_u64(reactant_raw),
                        )

                products = entry.get("products")
                if products:
                    for product_raw in products:
                        reaction_pairs += 1
                        p_hash = product_raw["hash"]
                        product_raw.pop("geo", None)
                        p_var_id: Optional[int] = None
                        if smiles_from_entry:
                            p_fp = _spec_fingerprint_u64(product_raw)
                            p_var_id = variant_id_by_fp_local.get(p_fp)
                            if p_var_id is None:
                                p_var_id = len(variant_specs_local)
                                variant_id_by_fp_local[p_fp] = p_var_id
                                variant_specs_local.append(product_raw)

                        if p_hash not in product_by_hash_local:
                            product_by_hash_local[p_hash] = product_raw
                            product_shape_sig_local[p_hash] = _spec_shape_sig(product_raw)
                        else:
                            prev_sig = product_shape_sig_local.get(p_hash)
                            cur_sig = _spec_shape_sig(product_raw)
                            if prev_sig is not None and prev_sig != cur_sig:
                                logger.warning(
                                    "Product record mismatch for hash=%s | prev=%s cur=%s | prev_fp=%016x cur_fp=%016x",
                                    p_hash,
                                    prev_sig,
                                    cur_sig,
                                    _spec_fingerprint_u64(product_by_hash_local[p_hash]),
                                    _spec_fingerprint_u64(product_raw),
                                )

                        if depth_filter and p_hash not in seen_species:
                            seen_species.add(p_hash)
                            next_frontier.add(p_hash)

                        # Hot loop: avoid function call + sorting for reaction_key.
                        key = (r_hash, p_hash) if r_hash <= p_hash else (p_hash, r_hash)
                        if key not in seen_reaction_keys_local:
                            seen_add(key)
                            # Keep first-seen orientation (A>>B vs B>>A) while deduplicating by unordered pair.
                            if smiles_from_entry:
                                if r_var_id is None or p_var_id is None:
                                    raise RuntimeError("Internal error: missing variant ids while smiles_from_entry.")
                                dedup_append((r_var_id, p_var_id))
                            else:
                                dedup_append((r_hash, p_hash))

                        if next_log_pairs is not None and reaction_pairs >= next_log_pairs:
                            now = time.time()
                            dt = now - last_log_t
                            dp = reaction_pairs - last_log_pairs
                            pairs_per_s = (dp / dt) if dt > 0 else 0.0
                            total_dt = now - read_start
                            total_pairs_per_s = (reaction_pairs / total_dt) if total_dt > 0 else 0.0
                            keep_pct = (100.0 * len(dedup_reactions) / reaction_pairs) if reaction_pairs else 0.0
                            rss_mb = _maxrss_mb()
                            rss_str = f"{rss_mb:.1f} MB" if rss_mb is not None else "n/a"
                            depth_str = f"{depth}" if depth_filter else "NA"
                            rem_str = "NA"
                            if depth_filter and depth >= 1 and frontier_remaining is not None:
                                rem_str = str(len(frontier_remaining))
                            unique_species = _union_len(reactant_by_hash_local, product_by_hash_local)
                            variants_str = ""
                            if smiles_from_entry:
                                variants_str = f" | variants={len(variant_specs_local) - 1}"
                            logger.info(
                                "Reading/Dedup | depth=%s | pairs=%d | unique_rxn=%d (%.1f%%) | unique_species=%d%s | rate=%.0f pairs/s (avg %.0f) | frontier_rem=%s | rss=%s",
                                depth_str,
                                reaction_pairs,
                                len(dedup_reactions),
                                keep_pct,
                                unique_species,
                                variants_str,
                                pairs_per_s,
                                total_pairs_per_s,
                                rem_str,
                                rss_str,
                            )
                            last_log_t = now
                            last_log_pairs = reaction_pairs
                            next_log_pairs += log_every_pairs

                if depth_filter and depth == 0 and depth0_remaining is not None and not depth0_remaining:
                    depth_elapsed = time.time() - depth_start_t
                    logger.info(
                        "Depth 0 done | reactants=%d | new_species=%d | elapsed %.2fs",
                        len(depth0_reactants),
                        len(next_frontier),
                        depth_elapsed,
                    )
                    if max_depth <= 0:
                        stop_reading = True
                        break
                    depth = 1
                    frontier_remaining = next_frontier
                    next_frontier = set()
                    frontier_size = len(frontier_remaining)
                    depth_start_t = time.time()
                    logger.info("Depth %d start | frontier=%d", depth, frontier_size)
                    if frontier_size == 0:
                        stop_reading = True
                        break

                if depth_filter and depth >= 1 and frontier_remaining is not None and not frontier_remaining:
                    depth_elapsed = time.time() - depth_start_t
                    logger.info(
                        "Depth %d done | reactants=%d | new_species=%d | elapsed %.2fs",
                        depth,
                        frontier_size,
                        len(next_frontier),
                        depth_elapsed,
                    )
                    if depth >= max_depth:
                        stop_reading = True
                        break
                    depth += 1
                    frontier_remaining = next_frontier
                    next_frontier = set()
                    frontier_size = len(frontier_remaining)
                    depth_start_t = time.time()
                    logger.info("Depth %d start | frontier=%d", depth, frontier_size)
                    if frontier_size == 0:
                        stop_reading = True
                        break

            if stop_reading:
                logger.info("Stopped reading early due to max_depth=%d", max_depth)

    read_elapsed = time.time() - run_start
    unique_species_total = _union_len(reactant_by_hash, product_by_hash)
    if smiles_from_entry:
        logger.info(
            "Read OK | reactants=%d | pairs=%d | unique_reactions=%d | unique_species=%d | variants=%d | elapsed %.2fs",
            reactant_records,
            reaction_pairs,
            len(dedup_reactions),
            unique_species_total,
            len(variant_specs) - 1,
            read_elapsed,
        )
    else:
        logger.info(
            "Read OK | reactants=%d | pairs=%d | unique_reactions=%d | unique_species=%d | elapsed %.2fs",
            reactant_records,
            reaction_pairs,
            len(dedup_reactions),
            unique_species_total,
            read_elapsed,
        )
    # Free the (potentially huge) seen-key set before RDKit work (stream mode only).
    if dedup_workers <= 1:
        seen_reaction_keys.clear()
        del seen_reaction_keys
        del seen_reaction_keys_local, seen_add

    # Build cleaned reactions and species set
    reactions_out: List[Dict] = []
    failed_smiles = 0

    total_reactions = len(dedup_reactions)
    if total_reactions == 0:
        logger.warning("No reactions found; writing empty output.")
    else:
        if max_workers <= 0:
            max_workers = os.cpu_count() or 1

        if smiles_from_entry:
            total_variants = len(variant_specs) - 1
            max_workers = max(1, min(max_workers, max(1, total_variants)))

            # Generate SMILES per *variant* (NOT per hash) to preserve atom-map nums per r->p entry.
            smiles_start = time.time()
            logger.info(
                "Generating variant SMILES (from entry records) | variants=%d | workers=%d | chunksize=%d",
                total_variants,
                max_workers,
                chunksize,
            )

            variant_smiles: Dict[int, Tuple[str, str]] = {}
            variant_tasks = (
                (var_id, spec["elements"], spec["bond_mats"])
                for var_id, spec in enumerate(variant_specs)
                if var_id and spec is not None
            )

            if max_workers > 1:
                ctx = mp.get_context("spawn")
                with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as executor:
                    v_iter: Iterable[Tuple[float, str, str, str]] = executor.map(
                        _species_smiles_task, variant_tasks, chunksize=max(1, chunksize)
                    )
                    if tqdm is not None:
                        v_iter = tqdm(
                            v_iter,
                            total=total_variants,
                            desc="Variant SMILES",
                            unit="variant",
                            disable=disable_tqdm,
                        )
                    for var_id, mapped, plain, error in v_iter:
                        if error:
                            failed_smiles += 1
                            if failed_smiles <= 10:
                                logger.warning("SMILES failed for variant %s: %s", var_id, error)
                        variant_smiles[int(var_id)] = (mapped, plain)
            else:
                v_iter2: Iterable[Tuple[float, str, str, str]] = (_species_smiles_task(task) for task in variant_tasks)
                if tqdm is not None:
                    v_iter2 = tqdm(
                        v_iter2,
                        total=total_variants,
                        desc="Variant SMILES",
                        unit="variant",
                        disable=disable_tqdm,
                    )
                for var_id, mapped, plain, error in v_iter2:
                    if error:
                        failed_smiles += 1
                        if failed_smiles <= 10:
                            logger.warning("SMILES failed for variant %s: %s", var_id, error)
                    variant_smiles[int(var_id)] = (mapped, plain)

            smiles_elapsed = time.time() - smiles_start
            logger.info("Variant SMILES done | elapsed %.2fs", smiles_elapsed)

            build_start = time.time()
            logger.info("Building reactions | reactions=%d", total_reactions)
            rxn_iter2 = dedup_reactions
            if tqdm is not None:
                rxn_iter2 = tqdm(
                    rxn_iter2,
                    total=total_reactions,
                    desc="Building reactions",
                    unit="reaction",
                    disable=disable_tqdm,
                )

            species_used: set[float] = set()
            for r_id, p_id in rxn_iter2:
                r_vid = int(r_id)
                p_vid = int(p_id)
                r_spec = variant_specs[r_vid]
                p_spec = variant_specs[p_vid]
                if r_spec is None or p_spec is None:
                    raise KeyError(f"Missing variant spec for ids {r_vid} -> {p_vid}")

                species_used.add(float(r_spec["hash"]))
                species_used.add(float(p_spec["hash"]))

                r_mapped, r_plain = variant_smiles.get(r_vid, ("", ""))
                p_mapped, p_plain = variant_smiles.get(p_vid, ("", ""))
                reactions_out.append(
                    {
                        "reactant": r_spec,
                        "product": p_spec,
                        "reaction_smiles": f"{r_plain}>>{p_plain}",
                        "reaction_smiles_mapped": f"{r_mapped}>>{p_mapped}",
                    }
                )

            logger.info("Unique species: %d", len(species_used))
            logger.info("Unique reactions: %d", len(reactions_out))
            if failed_smiles:
                logger.info("SMILES failures: %d (first 10 shown above)", failed_smiles)
            build_elapsed = time.time() - build_start
            logger.info("Build reactions done | elapsed %.2fs", build_elapsed)

        else:
            total_species_tasks = len(reactant_by_hash) + len(product_by_hash)
            max_workers = max(1, min(max_workers, max(1, total_species_tasks)))

            # Generate SMILES per-species role (NOT per-reaction) to avoid massive IPC and duplicated RDKit work.
            smiles_start = time.time()
            logger.info(
                "Generating species SMILES | reactant_species=%d | product_species=%d | workers=%d | chunksize=%d",
                len(reactant_by_hash),
                len(product_by_hash),
                max_workers,
                chunksize,
            )

            reactant_smiles: Dict[float, Tuple[str, str]] = {}
            product_smiles: Dict[float, Tuple[str, str]] = {}
            reactant_tasks = (
                (spec_hash, spec["elements"], spec["bond_mats"]) for spec_hash, spec in reactant_by_hash.items()
            )
            product_tasks = (
                (spec_hash, spec["elements"], spec["bond_mats"]) for spec_hash, spec in product_by_hash.items()
            )

            if max_workers > 1:
                # Use spawn to avoid forking a huge parent process (common immediate OOM-kill on HPC).
                ctx = mp.get_context("spawn")
                with ProcessPoolExecutor(max_workers=max_workers, mp_context=ctx) as executor:
                    rx_iter: Iterable[Tuple[float, str, str, str]] = executor.map(
                        _species_smiles_task, reactant_tasks, chunksize=max(1, chunksize)
                    )
                    if tqdm is not None:
                        rx_iter = tqdm(
                            rx_iter,
                            total=len(reactant_by_hash),
                            desc="Reactant SMILES",
                            unit="species",
                            disable=disable_tqdm,
                        )
                    for spec_hash, mapped, plain, error in rx_iter:
                        if error:
                            failed_smiles += 1
                            if failed_smiles <= 10:
                                logger.warning("SMILES failed for species %s: %s", spec_hash, error)
                        reactant_smiles[spec_hash] = (mapped, plain)

                    px_iter: Iterable[Tuple[float, str, str, str]] = executor.map(
                        _species_smiles_task, product_tasks, chunksize=max(1, chunksize)
                    )
                    if tqdm is not None:
                        px_iter = tqdm(
                            px_iter,
                            total=len(product_by_hash),
                            desc="Product SMILES",
                            unit="species",
                            disable=disable_tqdm,
                        )
                    for spec_hash, mapped, plain, error in px_iter:
                        if error:
                            failed_smiles += 1
                            if failed_smiles <= 10:
                                logger.warning("SMILES failed for species %s: %s", spec_hash, error)
                        product_smiles[spec_hash] = (mapped, plain)
            else:
                rx_iter2: Iterable[Tuple[float, str, str, str]] = (_species_smiles_task(task) for task in reactant_tasks)
                if tqdm is not None:
                    rx_iter2 = tqdm(
                        rx_iter2,
                        total=len(reactant_by_hash),
                        desc="Reactant SMILES",
                        unit="species",
                        disable=disable_tqdm,
                    )
                for spec_hash, mapped, plain, error in rx_iter2:
                    if error:
                        failed_smiles += 1
                        if failed_smiles <= 10:
                            logger.warning("SMILES failed for species %s: %s", spec_hash, error)
                    reactant_smiles[spec_hash] = (mapped, plain)

                px_iter2: Iterable[Tuple[float, str, str, str]] = (_species_smiles_task(task) for task in product_tasks)
                if tqdm is not None:
                    px_iter2 = tqdm(
                        px_iter2,
                        total=len(product_by_hash),
                        desc="Product SMILES",
                        unit="species",
                        disable=disable_tqdm,
                    )
                for spec_hash, mapped, plain, error in px_iter2:
                    if error:
                        failed_smiles += 1
                        if failed_smiles <= 10:
                            logger.warning("SMILES failed for species %s: %s", spec_hash, error)
                    product_smiles[spec_hash] = (mapped, plain)

            smiles_elapsed = time.time() - smiles_start
            logger.info("Species SMILES done | elapsed %.2fs", smiles_elapsed)

            build_start = time.time()
            logger.info("Building reactions | reactions=%d", total_reactions)
            rxn_iter: Iterable[Tuple[float, float]] = dedup_reactions
            if tqdm is not None:
                rxn_iter = tqdm(
                    rxn_iter,
                    total=total_reactions,
                    desc="Building reactions",
                    unit="reaction",
                    disable=disable_tqdm,
                )

            species_used: set[float] = set()
            missing_role_records = 0
            cross_role_fallbacks = 0
            for r_hash, p_hash in rxn_iter:
                r_spec = reactant_by_hash.get(r_hash)
                p_spec = product_by_hash.get(p_hash)
                if r_spec is None or p_spec is None:
                    missing_role_records += 1
                    if r_spec is None:
                        r_spec = product_by_hash.get(r_hash)
                        cross_role_fallbacks += 1 if r_spec is not None else 0
                    if p_spec is None:
                        p_spec = reactant_by_hash.get(p_hash)
                        cross_role_fallbacks += 1 if p_spec is not None else 0
                    if r_spec is None or p_spec is None:
                        raise KeyError(f"Missing species record for reaction {r_hash} -> {p_hash}")
                    if strict_role_records:
                        raise KeyError(
                            f"Strict role records enabled: missing role-specific record for reaction {r_hash} -> {p_hash}. "
                            "This can break atom-map correspondence; re-run without --strict-role-records to allow fallback."
                        )

                species_used.add(r_hash)
                species_used.add(p_hash)

                r_mapped, r_plain = reactant_smiles.get(r_hash, ("", ""))
                p_mapped, p_plain = product_smiles.get(p_hash, ("", ""))
                if (not r_mapped) and r_hash in product_smiles:
                    r_mapped, r_plain = product_smiles.get(r_hash, ("", ""))
                if (not p_mapped) and p_hash in reactant_smiles:
                    p_mapped, p_plain = reactant_smiles.get(p_hash, ("", ""))
                reactions_out.append(
                    {
                        "reactant": r_spec,
                        "product": p_spec,
                        "reaction_smiles": f"{r_plain}>>{p_plain}",
                        "reaction_smiles_mapped": f"{r_mapped}>>{p_mapped}",
                    }
                )

            logger.info("Unique species: %d", len(species_used))
            logger.info("Unique reactions: %d", len(reactions_out))
            if missing_role_records:
                logger.warning(
                    "Role-record misses while building reactions: %d (cross-role fallbacks=%d). "
                    "If atom map nums look wrong, try increasing --max-depth (to include missing reactants) or re-run with --strict-role-records to fail fast.",
                    missing_role_records,
                    cross_role_fallbacks,
                )
            if failed_smiles:
                logger.info("SMILES failures: %d (first 10 shown above)", failed_smiles)
            build_elapsed = time.time() - build_start
            logger.info("Build reactions done | elapsed %.2fs", build_elapsed)

    write_start = time.time()
    logger.info("Writing output pkl: %s", output_path)
    with open(output_path, "wb") as f:
        pickle.dump(reactions_out, f, protocol=pickle.HIGHEST_PROTOCOL)
    if verify_output_len:
        try:
            with open(output_path, "rb") as f2:
                obj = pickle.load(f2)
            if isinstance(obj, list):
                logger.info("Output verify | len=%d", len(obj))
            else:
                logger.warning("Output verify | loaded non-list object: %s", type(obj))
        except Exception as exc:
            logger.warning("Output verify failed: %s", exc)
    write_elapsed = time.time() - write_start
    total_elapsed = time.time() - run_start
    logger.info("Write OK | elapsed %.2fs", total_elapsed)
    logger.info(
        "Timing summary | read/dedup %.2fs | smiles %.2fs | build %.2fs | write %.2fs | total %.2fs",
        read_elapsed,
        smiles_elapsed,
        build_elapsed,
        write_elapsed,
        total_elapsed,
    )
    if total_elapsed > 0:
        read_frac = read_elapsed / total_elapsed
        smiles_frac = smiles_elapsed / total_elapsed
        if read_frac >= 0.45:
            logger.info(
                "Note: reading/dedup dominates runtime (%.0f%%). If still too slow, consider a bucketed multi-process dedup (plan 2).",
                read_frac * 100,
            )
        elif smiles_frac >= 0.45:
            logger.info(
                "Note: SMILES dominates runtime (%.0f%%). If you have spare CPUs, try higher --max-workers.",
                smiles_frac * 100,
            )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Post-process BNFM pkl to deduplicate and annotate reactions.")
    parser.add_argument("--input", required=True, help="Input pkl from bnfm_iterator.")
    parser.add_argument(
        "--output", default="bnfm_cleaned.pkl", help="Output pkl path for deduplicated reactions."
    )
    parser.add_argument(
        "--tqdm",
        action="store_true",
        help="Force enable tqdm progress bars (requires tqdm).",
    )
    parser.add_argument(
        "--no-tqdm",
        action="store_true",
        help="Disable tqdm progress bars.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        # default=0,
        default=8,
        help="Parallel workers for SMILES generation (0=auto). Use 1 to disable parallelism.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=50,
        help="ProcessPoolExecutor chunksize for SMILES generation (default: 50).",
    )
    parser.add_argument(
        "--log-every-pairs",
        type=int,
        default=1_000_000,
        help="Log progress every N reaction pairs during Reading/Dedup (0=disable).",
    )
    parser.add_argument(
        "--dedup-workers",
        type=int,
        default=1,
        help="Parallel workers for reaction-pair dedup during reading (1=disable).",
    )
    parser.add_argument(
        "--dedup-batch-pairs",
        type=int,
        default=50_000,
        help="Batch size (pairs) sent to each dedup worker (larger reduces IPC overhead).",
    )
    parser.add_argument(
        "--strict-role-records",
        action="store_true",
        help="Fail fast if a reaction needs cross-role fallback to find a species record (can affect atom map nums).",
    )
    parser.add_argument(
        "--verify-output-len",
        action="store_true",
        help="After writing output, reopen and log len(pickle.load(output)) (slow for huge outputs).",
    )
    parser.add_argument(
        "--smiles-from-entry",
        action="store_true",
        help="Generate reaction SMILES using the exact per-entry reactant/product records (variant-aware; preserves atom map nums across depths).",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=-1,
        help="Only parse reactions up to this reactant depth (0=seed reactants only; <0=disable). Requires seed SMILES.",
    )
    seed_group = parser.add_mutually_exclusive_group()
    seed_group.add_argument(
        "--seed-smiles",
        nargs="+",
        default=[],
        help="Depth-0 seed SMILES strings (whitespace-separated).",
    )
    seed_group.add_argument(
        "--seed-smiles-file",
        default="",
        help="Path to a text file containing depth-0 seed SMILES (whitespace-separated).",
    )
    seed_group.add_argument(
        "--seed-by-formula-file",
        default="",
        help="Path to smiles_by_formula.txt-style file (first token is formula, rest are SMILES). Use with --seed-formula or --seed-line.",
    )
    parser.add_argument(
        "--seed-formula",
        default="",
        help="Formula key to select a line from --seed-by-formula-file.",
    )
    parser.add_argument(
        "--seed-line",
        type=int,
        default=0,
        help="1-based line number to select from --seed-by-formula-file (alternative to --seed-formula).",
    )
    args = parser.parse_args()
    if args.max_depth is not None and args.max_depth >= 0:
        if args.seed_by_formula_file and not (args.seed_formula or args.seed_line):
            parser.error("--seed-by-formula-file requires --seed-formula or --seed-line when --max-depth >= 0.")
        if not (args.seed_smiles or args.seed_smiles_file or args.seed_by_formula_file):
            parser.error(
                "--max-depth >= 0 requires seed SMILES: provide --seed-smiles, --seed-smiles-file, or --seed-by-formula-file (+ --seed-formula/--seed-line)."
            )
    return args


def main() -> None:
    args = _parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    use_tqdm = (tqdm is not None) and (not args.no_tqdm)
    if args.tqdm and (tqdm is None):
        logger.warning("--tqdm requested but tqdm is not installed; continuing without progress bars.")
    if args.tqdm:
        use_tqdm = True

    depth0_seed_hashes: Optional[set[float]] = None
    if args.max_depth is not None and args.max_depth >= 0:
        if args.seed_smiles:
            seed_smiles = args.seed_smiles
        elif args.seed_smiles_file:
            seed_smiles = _load_seed_smiles_file(args.seed_smiles_file)
        elif args.seed_by_formula_file:
            if args.seed_line:
                formula, seed_smiles = _load_seed_smiles_by_line(args.seed_by_formula_file, args.seed_line)
                if args.seed_formula and args.seed_formula != formula:
                    logger.warning(
                        "seed_formula=%s does not match formula on line %d (%s); using line selection",
                        args.seed_formula,
                        args.seed_line,
                        formula,
                    )
            else:
                seed_smiles = _load_seed_smiles_by_formula(args.seed_by_formula_file, args.seed_formula)
        else:
            seed_smiles = []
        depth0_seed_hashes = _hashes_from_smiles(seed_smiles)
        if not depth0_seed_hashes:
            raise ValueError("No valid seed hashes computed from the provided seed SMILES.")

    process(
        args.input,
        args.output,
        use_tqdm=use_tqdm,
        max_workers=args.max_workers,
        chunksize=args.chunksize,
        log_every_pairs=args.log_every_pairs,
        dedup_workers=args.dedup_workers,
        dedup_batch_pairs=args.dedup_batch_pairs,
        max_depth=args.max_depth,
        depth0_seed_hashes=depth0_seed_hashes,
        strict_role_records=args.strict_role_records,
        verify_output_len=args.verify_output_len,
        smiles_from_entry=args.smiles_from_entry,
    )


if __name__ == "__main__":
    main()
