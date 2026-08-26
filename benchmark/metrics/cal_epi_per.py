#!/usr/bin/env python3
"""
Compare epitope retention between two HCDR3-antigen complexes.

Given two PDB structures of HCDR3-antigen complexes, this script:
  1. Identifies epitope residues in complex 1 (antigen surface residues whose
     C-beta distance to any HCDR3 residue is below the cutoff).
  2. Maps those epitope residues onto complex 2 (optionally after sequence +
     structure alignment, useful when complex 2 is a truncated variant).
  3. Reports the fraction of epitope residues retained in complex 2.

Requirements:
    pip install biopython numpy
"""

import argparse
import csv
import os
import sys
import numpy as np

# Fix encoding for Windows GBK terminals
if sys.stdout.encoding and sys.stdout.encoding.upper() in ("GBK", "CP936"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# ── Biopython imports ──────────────────────────────────────────────────────
try:
    from Bio.PDB import PDBParser, Superimposer
    from Bio.PDB.Polypeptide import protein_letters_3to1
    from Bio.Align import PairwiseAligner, substitution_matrices
except ImportError:
    sys.exit(
        "This script requires Biopython.  Install it with:\n"
        "    pip install biopython"
    )


# ═══════════════════════════════════════════════════════════════════════════
#  Helper utilities
# ═══════════════════════════════════════════════════════════════════════════

def _get_cb_atom(residue):
    """Return the CB atom of *residue*, falling back to CA for glycine."""
    if "CB" in residue:
        return residue["CB"]
    # Glycine (and other edge cases) → use CA
    if "CA" in residue:
        return residue["CA"]
    return None


def _get_residue_sequence(residues):
    """Convert a list of Bio.PDB Residue objects to a single-letter sequence string.

    Non-standard residues are mapped to 'X'.
    """
    seq_parts = []
    for res in residues:
        resname = res.get_resname().strip()
        seq_parts.append(protein_letters_3to1.get(resname, "X"))
    return "".join(seq_parts)


def _get_chain_residues(structure, chain_id):
    """Return ordered (hetflag=False, standard) residues for *chain_id* in the
    first model of *structure*."""
    model = structure[0]
    if chain_id not in model:
        available = ", ".join(model.child_dict.keys())
        raise ValueError(
            f"Chain '{chain_id}' not found in {structure.get_id()}. "
            f"Available chains: {available}"
        )
    chain = model[chain_id]
    # Filter to standard residues (HETATM excluded, water excluded)
    residues = [
        res
        for res in chain.get_residues()
        if res.get_id()[0] == " "  # standard residue flag
    ]
    return residues


# ═══════════════════════════════════════════════════════════════════════════
#  Epitope detection
# ═══════════════════════════════════════════════════════════════════════════

def compute_epitope_residues(hcdr3_residues, antigen_residues, cutoff):
    """Identify epitope residues and their distances to HCDR3 in complex 1.

    Parameters
    ----------
    hcdr3_residues : list of Bio.PDB.Residue
    antigen_residues : list of Bio.PDB.Residue
    cutoff : float
        Distance threshold in Angstrom.

    Returns
    -------
    epitope_residues : list of Bio.PDB.Residue
    distances : list of float
        Minimum C-beta distance from each epitope residue to any HCDR3 residue.
    """
    # Pre-compute CB coordinate for every HCDR3 residue
    hcdr3_coords = []
    for res in hcdr3_residues:
        cb = _get_cb_atom(res)
        if cb is not None:
            hcdr3_coords.append(cb.get_coord())
    if not hcdr3_coords:
        return [], []

    hcdr3_coords = np.array(hcdr3_coords)  # (N, 3)

    epitope = []
    distances = []
    for ag_res in antigen_residues:
        cb = _get_cb_atom(ag_res)
        if cb is None:
            continue
        coord = cb.get_coord()  # (3,)
        dists = np.linalg.norm(hcdr3_coords - coord, axis=1)
        min_d = np.min(dists)
        if min_d < cutoff:
            epitope.append(ag_res)
            distances.append(float(min_d))

    return epitope, distances


# ═══════════════════════════════════════════════════════════════════════════
#  Sequence alignment
# ═══════════════════════════════════════════════════════════════════════════

def align_sequences(seq1, seq2):
    """Perform local pairwise alignment of two protein sequences.

    Returns
    -------
    aligned1 : str
        Aligned sequence 1 (with gaps).
    aligned2 : str
        Aligned sequence 2 (with gaps).
    mapping : dict
        Maps 0-based position in seq1 → 0-based position in seq2 (only for
        aligned columns where both residues are not gaps).
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    alignments = aligner.align(seq1, seq2)
    if not alignments:
        raise RuntimeError("Sequence alignment returned no results.")

    best = alignments[0]
    # Direct indexing: best[i] returns the i-th aligned sequence (with gaps)
    aligned1 = str(best[0])
    aligned2 = str(best[1])

    # Build position mapping: position in original seq1 → position in original seq2
    mapping = {}
    pos1 = 0
    pos2 = 0
    for a, b in zip(aligned1, aligned2):
        if a != "-" and b != "-":
            mapping[pos1] = pos2
        if a != "-":
            pos1 += 1
        if b != "-":
            pos2 += 1

    return aligned1, aligned2, mapping


# ═══════════════════════════════════════════════════════════════════════════
#  Structure superposition
# ═══════════════════════════════════════════════════════════════════════════

def superpose_structures(fixed_residues, moving_residues, moving_structure):
    """Superpose *moving_residues* onto *fixed_residues* using C-alpha atoms,
    then apply the transformation to every atom in *moving_structure*.

    Parameters
    ----------
    fixed_residues : list of Bio.PDB.Residue
        Reference antigen residues (complex 1).
    moving_residues : list of Bio.PDB.Residue
        Antigen residues to be moved (complex 2).  Must be 1:1 with fixed_residues
        after sequence alignment.
    moving_structure : Bio.PDB.Structure.Structure
        The entire second complex structure — all its atoms will be transformed.

    Returns
    -------
    float
        RMSD of the superposition (Å).
    """
    # Collect CA atoms from matched residue pairs
    fixed_atoms = []
    moving_atoms = []
    for fr, mr in zip(fixed_residues, moving_residues):
        if "CA" in fr and "CA" in mr:
            fixed_atoms.append(fr["CA"])
            moving_atoms.append(mr["CA"])

    if len(fixed_atoms) < 3:
        raise RuntimeError(
            f"Need at least 3 matched C-alpha pairs for superposition; "
            f"found {len(fixed_atoms)}."
        )

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)
    rot, tran = sup.rotran

    # Apply transformation to every atom in the moving structure
    for model in moving_structure:
        for chain in model:
            for atom in chain.get_atoms():
                atom.transform(rot, tran)

    return sup.rms


# ═══════════════════════════════════════════════════════════════════════════
#  Retention calculation
# ═══════════════════════════════════════════════════════════════════════════

def compute_retention(
    epitope_residues,
    epitope_distances1,
    ag1_residues,
    ag2_residues,
    hcdr3_2_residues,
    cutoff,
    residue_mapping=None,
):
    """Compute the fraction of epitope residues retained in the second complex.

    A residue is "retained" if:
      - It can be mapped to a residue in the second antigen (via
        *residue_mapping* if provided, otherwise by matching residue number).
      - Its C-beta distance to the second HCDR3 is still below *cutoff*.

    Parameters
    ----------
    epitope_residues : list of Bio.PDB.Residue
        Epitope residues from complex 1.
    epitope_distances1 : list of float
        Min CB distance from each epitope residue to HCDR3 in complex 1.
    ag1_residues : list of Bio.PDB.Residue
        All antigen residues from complex 1 (used to resolve index in alignment mapping).
    ag2_residues : list of Bio.PDB.Residue
        All antigen residues from complex 2.
    hcdr3_2_residues : list of Bio.PDB.Residue
        HCDR3 residues from complex 2.
    cutoff : float
        Distance threshold.
    residue_mapping : dict | None
        If provided, maps 0-based index in ag1 → 0-based index in ag2.
        Otherwise, residues are matched by residue number (res.get_id()[1]).

    Returns
    -------
    fraction_retained : float
    retained_details : list of dict
        Per-epitope-residue information.  Keys: epitope_index, resname, resnum,
        dist1, dist2, retained, match_method, reason.
    """
    # Build lookups: (resnum, insertion_code) → index / (index, residue)
    ag1_by_key = {}
    for i, res in enumerate(ag1_residues):
        key = (res.get_id()[1], res.get_id()[2])
        ag1_by_key[key] = i

    ag2_by_key = {}
    for i, res in enumerate(ag2_residues):
        key = (res.get_id()[1], res.get_id()[2])
        ag2_by_key[key] = (i, res)

    # Pre-compute HCDR3 CB coords for complex 2
    hcdr3_2_coords = []
    for res in hcdr3_2_residues:
        cb = _get_cb_atom(res)
        if cb is not None:
            hcdr3_2_coords.append(cb.get_coord())
    hcdr3_2_coords = np.array(hcdr3_2_coords) if hcdr3_2_coords else np.empty((0, 3))

    retained_count = 0
    details = []

    for epi_idx, epi_res in enumerate(epitope_residues):
        epi_resnum = epi_res.get_id()[1]
        epi_resname = epi_res.get_resname().strip()

        # Determine the corresponding residue in complex 2
        matched_res = None
        match_method = ""

        # Resolve this epitope residue's index within ag1
        epi_key = (epi_res.get_id()[1], epi_res.get_id()[2])
        ag1_idx = ag1_by_key.get(epi_key)
        if ag1_idx is None:
            match_method = f"residue {epi_resnum} not found in ag1 index"
        elif residue_mapping is not None:
            # Use sequence-alignment mapping (ag1 index → ag2 index)
            if ag1_idx in residue_mapping:
                ag2_idx = residue_mapping[ag1_idx]
                if ag2_idx < len(ag2_residues):
                    matched_res = ag2_residues[ag2_idx]
                    match_method = f"alignment (ag1[{ag1_idx}]→ag2[{ag2_idx}])"
            else:
                match_method = f"no alignment mapping for ag1 index {ag1_idx}"
        else:
            # Match by (resnum, insertion_code) key
            if epi_key in ag2_by_key:
                _, matched_res = ag2_by_key[epi_key]
                match_method = f"residue key {epi_key}"

        # Check retention
        retained = False
        reason = ""
        dist2 = None

        if matched_res is None:
            reason = "not present in complex 2 antigen"
        else:
            cb2 = _get_cb_atom(matched_res)
            if cb2 is None:
                reason = "no CB/CA atom in complex 2"
            elif len(hcdr3_2_coords) == 0:
                reason = "no HCDR3 atoms in complex 2"
            else:
                dist2 = float(np.min(
                    np.linalg.norm(hcdr3_2_coords - cb2.get_coord(), axis=1)
                ))
                if dist2 < cutoff:
                    retained = True
                    retained_count += 1
                    reason = f"distance {dist2:.2f} Å < {cutoff:.1f} Å"
                else:
                    reason = f"distance {dist2:.2f} Å ≥ {cutoff:.1f} Å"

        dist1 = epitope_distances1[epi_idx]

        details.append(
            {
                "epitope_index": epi_idx,
                "resname": epi_resname,
                "resnum": epi_resnum,
                "dist1": dist1,
                "dist2": dist2,
                "retained": retained,
                "match_method": match_method,
                "reason": reason,
            }
        )

    fraction = retained_count / len(epitope_residues) if epitope_residues else 0.0
    return fraction, details


# ═══════════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="Compare epitope retention between two HCDR3-antigen complexes."
    )
    # Positional arguments
    parser.add_argument("pdb1", help="Path to first HCDR3-antigen complex PDB file")
    parser.add_argument("pdb2", help="Path to second HCDR3-antigen complex PDB file")
    parser.add_argument("hcdr3_chain_1", help="Chain ID of HCDR3 in first complex")
    parser.add_argument("antigen_chain_1", help="Chain ID of antigen in first complex")
    parser.add_argument("hcdr3_chain_2", help="Chain ID of HCDR3 in second complex")
    parser.add_argument("antigen_chain_2", help="Chain ID of antigen in second complex")

    # Optional named arguments
    parser.add_argument(
        "-c", "--epitope-cutoff",
        type=float,
        default=8.0,
        help="Distance cutoff (Å) for epitope definition. "
             "An antigen residue is epitope if its C-beta distance "
             "to any HCDR3 residue is below this threshold (default: 8.0).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print detailed per-residue retention information.",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Path to output CSV file. If the file exists, a new row is appended.",
    )
    parser.add_argument(
        "-a", "--align",
        action="store_true",
        help="After sequence-based renumbering, also perform structural "
             "superposition (rotate + translate complex 2 onto complex 1). "
             "Use this when the two structures may be in different "
             "coordinate frames or one is a truncated variant of the other. "
             "Without this flag, only sequence alignment (renumbering) is done.",
    )

    args = parser.parse_args()

    # ── Basic validation ───────────────────────────────────────────────
    for pdb_path in (args.pdb1, args.pdb2):
        if not os.path.isfile(pdb_path):
            sys.exit(f"PDB file not found: {pdb_path}")

    # ── Parse structures ───────────────────────────────────────────────
    parser_pdb = PDBParser(QUIET=True)
    struct1 = parser_pdb.get_structure("complex1", args.pdb1)
    struct2 = parser_pdb.get_structure("complex2", args.pdb2)

    # ── Extract residues ───────────────────────────────────────────────
    try:
        hcdr3_1 = _get_chain_residues(struct1, args.hcdr3_chain_1)
        ag1 = _get_chain_residues(struct1, args.antigen_chain_1)
        hcdr3_2 = _get_chain_residues(struct2, args.hcdr3_chain_2)
        ag2 = _get_chain_residues(struct2, args.antigen_chain_2)
    except ValueError as e:
        sys.exit(str(e))

    print(f"Complex 1 — HCDR3 chain {args.hcdr3_chain_1}: {len(hcdr3_1)} residues, "
          f"antigen chain {args.antigen_chain_1}: {len(ag1)} residues")
    print(f"Complex 2 — HCDR3 chain {args.hcdr3_chain_2}: {len(hcdr3_2)} residues, "
          f"antigen chain {args.antigen_chain_2}: {len(ag2)} residues")

    # ── Epitope in complex 1 ───────────────────────────────────────────
    epitope, epitope_dist1 = compute_epitope_residues(hcdr3_1, ag1, args.epitope_cutoff)
    print(f"\nEpitope residues in complex 1 (cutoff < {args.epitope_cutoff} Å): "
          f"{len(epitope)} / {len(ag1)}")

    if not epitope:
        print("WARNING: No epitope residues found in complex 1. "
              "Retention fraction is undefined (set to 0).")
        fraction = 0.0
        details = []
    else:
        # ── Sequence alignment (renumbering) — always performed ────────
        # Even when the antigen is the same protein, residue/atom numbering
        # may differ between files, so we always align to map residues.
        print("\n── Performing sequence alignment (renumbering) ──")
        seq1 = _get_residue_sequence(ag1)
        seq2 = _get_residue_sequence(ag2)
        print(f"Antigen 1 sequence length: {len(seq1)}")
        print(f"Antigen 2 sequence length: {len(seq2)}")

        aln1, aln2, residue_mapping = align_sequences(seq1, seq2)
        n_matched = len(residue_mapping)
        print(f"Aligned positions matched: {n_matched} / {len(seq1)}")

        # ── Optional structural superposition ──────────────────────────
        if args.align:
            print("\n── Performing structural superposition ──")
            fixed_list = []
            moving_list = []
            for i1, i2 in sorted(residue_mapping.items()):
                fixed_list.append(ag1[i1])
                moving_list.append(ag2[i2])

            rmsd = superpose_structures(fixed_list, moving_list, struct2)
            print(f"Superposition RMSD: {rmsd:.3f} Å "
                  f"(over {len(fixed_list)} C-alpha pairs)")

        # ── Compute retention ──────────────────────────────────────────
        fraction, details = compute_retention(
            epitope_residues=epitope,
            epitope_distances1=epitope_dist1,
            ag1_residues=ag1,
            ag2_residues=ag2,
            hcdr3_2_residues=hcdr3_2,
            cutoff=args.epitope_cutoff,
            residue_mapping=residue_mapping,
        )

    # ── Result summary ─────────────────────────────────────────────────
    retained_count = sum(1 for d in details if d["retained"])
    print(f"\n{'='*60}")
    print(f"Epitope retention: {retained_count} / {len(epitope)} "
          f"({fraction:.2%})")
    print(f"{'='*60}")

    # ── Verbose output ─────────────────────────────────────────────────
    if args.verbose and details:
        print(f"\n── Per-residue details ──")
        header = (f"{'Idx':>4s}  {'Res':>5s}  {'Num':>5s}  "
                  f"{'Dist1':>7s}  {'Dist2':>7s}  {'Retained':>8s}  "
                  f"{'Match method':<30s}  {'Reason'}")
        print(header)
        print("-" * len(header))
        for d in details:
            d1 = f"{d['dist1']:.2f}" if d['dist1'] is not None else "N/A"
            d2 = f"{d['dist2']:.2f}" if d['dist2'] is not None else "N/A"
            print(
                f"{d['epitope_index']:4d}  "
                f"{d['resname']:>5s}  "
                f"{d['resnum']:5d}  "
                f"{d1:>7s}  "
                f"{d2:>7s}  "
                f"{'YES' if d['retained'] else 'NO':>8s}  "
                f"{d['match_method']:<30s}  "
                f"{d['reason']}"
            )

    # ── CSV output ─────────────────────────────────────────────────────
    if args.output:
        file_exists = os.path.isfile(args.output)
        with open(args.output, "a", newline="") as fh:
            writer = csv.writer(fh)
            if not file_exists:
                writer.writerow(
                    ["PDB1", "PDB2", "epitope_residue",
                     "distance_complex1", "distance_complex2", "retained",
                     "retention_fraction"]
                )
            for d in details:
                label = f"{d['resname']}{d['resnum']}"
                d1 = f"{d['dist1']:.3f}" if d['dist1'] is not None else "N/A"
                d2 = f"{d['dist2']:.3f}" if d['dist2'] is not None else "N/A"
                writer.writerow(
                    [args.pdb1, args.pdb2, label, d1, d2,
                     "YES" if d["retained"] else "NO",
                     f"{fraction:.4f}"]
                )
        print(f"\nPer-residue results written to: {args.output}")


if __name__ == "__main__":
    main()

