import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

from enumerate_reaction_mapping import generate_atom_mapped_reactions


def _compute_mappings(rsmi):
    rsmi = (rsmi or "").strip()
    if not rsmi:
        return []
    return generate_atom_mapped_reactions(rxn_smiles=rsmi)


def _default_checkpoint_path(output_csv):
    return f"{output_csv}.ckpt.json"


def _input_signature(path):
    stat = os.stat(path)
    return {"size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _load_checkpoint(path):
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except FileNotFoundError:
        return None
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid checkpoint JSON: {path}") from e


def _save_checkpoint(path, state):
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    tmp_path = f"{path}.tmp"
    with open(tmp_path, "w", encoding="utf-8") as handle:
        json.dump(state, handle, ensure_ascii=False, sort_keys=True, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp_path, path)


def main():
    parser = argparse.ArgumentParser(
        description="Augment reaction SMILES in a CSV using atom-mapped enumeration."
    )
    parser.add_argument("input_csv", help="Input CSV path.")
    parser.add_argument("output_csv", help="Output CSV path.")
    parser.add_argument(
        "--origin-idx-col",
        default="origin_idx",
        help="Column name for origin_idx (default: origin_idx).",
    )
    parser.add_argument(
        "--origin-rsmi-col",
        default="origin_rsmi",
        help="Column name for origin_rsmi (default: origin_rsmi).",
    )
    parser.add_argument(
        "--aug-idx-col",
        default="aug_idx",
        help="Column name for augmented index (default: aug_idx).",
    )
    parser.add_argument(
        "--keep-empty",
        action="store_true",
        help="Keep rows with no mappings (aug_idx=0).",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print per-row progress to stderr.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of worker processes; 0 uses all CPUs.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from an existing output+checkpoint instead of overwriting.",
    )
    parser.add_argument(
        "--checkpoint",
        default=None,
        help="Checkpoint JSON path (default: <output_csv>.ckpt.json).",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=1,
        help="Write checkpoint every N input rows (default: 1).",
    )
    parser.add_argument(
        "--flush-every",
        type=int,
        default=1,
        help="Flush output every N input rows (default: 1).",
    )
    parser.add_argument(
        "--fsync",
        action="store_true",
        help="Also fsync output on each flush (slower, safer).",
    )
    args = parser.parse_args()

    if args.checkpoint_every < 1:
        raise ValueError("--checkpoint-every must be >= 1.")
    if args.flush_every < 1:
        raise ValueError("--flush-every must be >= 1.")

    checkpoint_path = args.checkpoint or _default_checkpoint_path(args.output_csv)
    input_sig = _input_signature(args.input_csv)

    checkpoint = _load_checkpoint(checkpoint_path)
    completed_rows = 0

    if args.resume:
        if checkpoint is None:
            raise ValueError(
                f"--resume requires an existing checkpoint: {checkpoint_path}"
            )
        if checkpoint.get("input_csv") != args.input_csv:
            raise ValueError(
                f"Checkpoint input_csv mismatch: {checkpoint.get('input_csv')} != {args.input_csv}"
            )
        if checkpoint.get("output_csv") != args.output_csv:
            raise ValueError(
                f"Checkpoint output_csv mismatch: {checkpoint.get('output_csv')} != {args.output_csv}"
            )
        if checkpoint.get("origin_idx_col") != args.origin_idx_col:
            raise ValueError("--resume with different --origin-idx-col is not supported.")
        if checkpoint.get("origin_rsmi_col") != args.origin_rsmi_col:
            raise ValueError(
                "--resume with different --origin-rsmi-col is not supported."
            )
        if checkpoint.get("aug_idx_col") != args.aug_idx_col:
            raise ValueError("--resume with different --aug-idx-col is not supported.")
        if checkpoint.get("keep_empty") != args.keep_empty:
            raise ValueError("--resume with different --keep-empty is not supported.")
        if checkpoint.get("input_signature") != input_sig:
            raise ValueError(
                "Input CSV has changed since the checkpoint was created; "
                "refuse to resume to avoid corrupt output."
            )
        if not os.path.exists(args.output_csv):
            raise ValueError(
                f"--resume requires an existing output file: {args.output_csv}"
            )
        if os.path.getsize(args.output_csv) == 0:
            raise ValueError(
                f"--resume output file is empty/corrupt: {args.output_csv}"
            )
        completed_rows = int(checkpoint.get("completed_rows", 0) or 0)
        if args.debug:
            print(
                f"[debug] resume from row={completed_rows + 1} (checkpoint={checkpoint_path})",
                file=sys.stderr,
            )
    else:
        checkpoint = {
            "version": 1,
            "input_csv": args.input_csv,
            "output_csv": args.output_csv,
            "origin_idx_col": args.origin_idx_col,
            "origin_rsmi_col": args.origin_rsmi_col,
            "aug_idx_col": args.aug_idx_col,
            "keep_empty": bool(args.keep_empty),
            "input_signature": input_sig,
            "created_at": time.time(),
            "completed_rows": 0,
            "updated_at": time.time(),
        }
        _save_checkpoint(checkpoint_path, checkpoint)

    with open(args.input_csv, "r", encoding="utf-8", newline="") as infile:
        reader = csv.DictReader(infile)
        if reader.fieldnames is None:
            raise ValueError("Input CSV has no header.")
        if args.origin_idx_col not in reader.fieldnames:
            raise ValueError(f"Missing column: {args.origin_idx_col}")
        if args.origin_rsmi_col not in reader.fieldnames:
            raise ValueError(f"Missing column: {args.origin_rsmi_col}")

        fieldnames = list(reader.fieldnames)
        if args.aug_idx_col not in fieldnames:
            fieldnames.append(args.aug_idx_col)

        output_mode = "a" if args.resume else "w"
        with open(args.output_csv, output_mode, encoding="utf-8", newline="") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            if not args.resume:
                writer.writeheader()

            rows_written_since_flush = 0
            rows_written_since_ckpt = 0
            last_completed_row = completed_rows

            def mark_row_done(row_num):
                nonlocal rows_written_since_flush, rows_written_since_ckpt, last_completed_row

                last_completed_row = row_num
                rows_written_since_flush += 1
                rows_written_since_ckpt += 1

                # Never allow the checkpoint to get ahead of what has been flushed to disk.
                if (
                    rows_written_since_flush >= args.flush_every
                    or rows_written_since_ckpt >= args.checkpoint_every
                ):
                    outfile.flush()
                    if args.fsync:
                        os.fsync(outfile.fileno())
                    rows_written_since_flush = 0

                if rows_written_since_ckpt >= args.checkpoint_every:
                    checkpoint["completed_rows"] = row_num
                    checkpoint["updated_at"] = time.time()
                    _save_checkpoint(checkpoint_path, checkpoint)
                    rows_written_since_ckpt = 0

            def write_outputs(row_num, row, rsmi, mappings):
                if not rsmi:
                    if args.keep_empty:
                        out_row = dict(row)
                        out_row[args.aug_idx_col] = "0"
                        writer.writerow(out_row)
                    if args.debug:
                        print(f"[debug] row={row_num} empty rsmi", file=sys.stderr)
                    return

                if args.debug:
                    origin_idx = row.get(args.origin_idx_col, "")
                    print(
                        f"[debug] row={row_num} origin_idx={origin_idx} mappings={len(mappings)}",
                        file=sys.stderr,
                    )

                if not mappings:
                    if args.keep_empty:
                        out_row = dict(row)
                        out_row[args.origin_rsmi_col] = rsmi
                        out_row[args.aug_idx_col] = "0"
                        writer.writerow(out_row)
                    return

                for idx, mapped in enumerate(mappings, start=1):
                    out_row = dict(row)
                    out_row[args.origin_rsmi_col] = mapped
                    out_row[args.aug_idx_col] = str(idx)
                    writer.writerow(out_row)

            jobs = args.jobs
            if jobs == 0:
                jobs = os.cpu_count() or 1
            if jobs < 1:
                raise ValueError("--jobs must be >= 1 (or 0 for all CPUs).")

            if jobs == 1:
                for row_num, row in enumerate(reader, start=1):
                    if row_num <= completed_rows:
                        continue
                    rsmi = (row.get(args.origin_rsmi_col) or "").strip()
                    mappings = _compute_mappings(rsmi)
                    write_outputs(row_num, row, rsmi, mappings)
                    mark_row_done(row_num)
                # Ensure final flush + checkpoint.
                outfile.flush()
                if args.fsync:
                    os.fsync(outfile.fileno())
                checkpoint["completed_rows"] = last_completed_row
                checkpoint["updated_at"] = time.time()
                _save_checkpoint(checkpoint_path, checkpoint)
                return

            pending = {}
            buffer = {}
            next_row_to_write = completed_rows + 1
            max_inflight = max(1, jobs * 4)

            with ProcessPoolExecutor(max_workers=jobs) as executor:
                for row_num, row in enumerate(reader, start=1):
                    if row_num <= completed_rows:
                        continue
                    rsmi = (row.get(args.origin_rsmi_col) or "").strip()
                    if not rsmi:
                        buffer[row_num] = (row, rsmi, [])
                    else:
                        fut = executor.submit(_compute_mappings, rsmi)
                        pending[fut] = (row_num, row, rsmi)

                    while pending and len(pending) >= max_inflight:
                        done, _ = wait(pending, return_when=FIRST_COMPLETED)
                        for fut in done:
                            rnum, r, rsm = pending.pop(fut)
                            try:
                                mappings = fut.result()
                            except Exception as e:
                                origin_idx = r.get(args.origin_idx_col, "")
                                raise RuntimeError(
                                    f"Failed mapping at row={rnum} origin_idx={origin_idx}"
                                ) from e
                            buffer[rnum] = (r, rsm, mappings)

                        while next_row_to_write in buffer:
                            r, rsm, mappings = buffer.pop(next_row_to_write)
                            write_outputs(next_row_to_write, r, rsm, mappings)
                            mark_row_done(next_row_to_write)
                            next_row_to_write += 1

                while pending:
                    done, _ = wait(pending, return_when=FIRST_COMPLETED)
                    for fut in done:
                        rnum, r, rsm = pending.pop(fut)
                        try:
                            mappings = fut.result()
                        except Exception as e:
                            origin_idx = r.get(args.origin_idx_col, "")
                            raise RuntimeError(
                                f"Failed mapping at row={rnum} origin_idx={origin_idx}"
                            ) from e
                        buffer[rnum] = (r, rsm, mappings)

                    while next_row_to_write in buffer:
                        r, rsm, mappings = buffer.pop(next_row_to_write)
                        write_outputs(next_row_to_write, r, rsm, mappings)
                        mark_row_done(next_row_to_write)
                        next_row_to_write += 1

            while next_row_to_write in buffer:
                r, rsm, mappings = buffer.pop(next_row_to_write)
                write_outputs(next_row_to_write, r, rsm, mappings)
                mark_row_done(next_row_to_write)
                next_row_to_write += 1

            # Ensure final flush + checkpoint.
            outfile.flush()
            if args.fsync:
                os.fsync(outfile.fileno())
            checkpoint["completed_rows"] = last_completed_row
            checkpoint["updated_at"] = time.time()
            _save_checkpoint(checkpoint_path, checkpoint)


if __name__ == "__main__":
    main()
