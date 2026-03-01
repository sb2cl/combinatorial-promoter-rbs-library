#!/usr/bin/env python3
"""Validate that the paper table matches the GenBank library.

Checks:
- Every table row maps to an existing GenBank file (by filename column if present, else by heuristic)
- Promoter and RBS labels in the GenBank FEATURES match the table
- Promoter feature sequence is 35 bp (J231xx set)
- RBS is annotated as regulatory with regulatory_class=ribosome_binding_site
- CDS translations have no internal stops (basic check)

Exit code:
- 0 if all checks pass
- 1 if any check fails
"""

import os, re, glob, argparse
import pandas as pd
from Bio import SeqIO

def sanitize_locus(text: str) -> str:
    lines = text.splitlines(True)
    if not lines or not lines[0].startswith("LOCUS"):
        return text
    line = lines[0].rstrip("\n")
    m = re.search(r"LOCUS\s+(\S+)\s+(\d+)\s+bp", line)
    if not m:
        return text
    name = m.group(1)
    length = int(m.group(2))
    mol = "ds-DNA" if "ds-DNA" in line else "DNA"
    topo = "circular" if "circular" in line else ("linear" if "linear" in line else "circular")
    div = "SYN" if "SYN" in line else "UNK"
    dm = re.search(r"(\d{2}-[A-Z]{3}-\d{4})", line)
    dt = dm.group(1) if dm else "01-JAN-1980"
    lines[0] = f"LOCUS       {name:<16}{length:>11} bp    {mol:<6} {topo:<8} {div:<3} {dt}\n"
    return "".join(lines)

def parse_genbank_relaxed(path):
    raw = open(path, "r", encoding="utf-8", errors="replace").read()
    fixed = sanitize_locus(raw)
    import tempfile
    with tempfile.NamedTemporaryFile("w+", delete=False, suffix=".gb") as tf:
        tf.write(fixed)
        tmp = tf.name
    try:
        rec = next(SeqIO.parse(tmp, "genbank"))
    finally:
        os.unlink(tmp)
    return rec

def get_qual(feat, key):
    vals = feat.qualifiers.get(key)
    return vals[0] if vals else None

def norm(x): 
    return str(x).strip().replace("-", "_")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="metadata/library_table_from_paper.csv")
    ap.add_argument("--genbank-dir", default="genbank")
    args = ap.parse_args()

    df = pd.read_csv(args.table)
    df["Promoter_norm"] = df["Promoter"].map(norm)
    df["RBS_norm"] = df["RBS"].map(norm)

    gb_paths = glob.glob(os.path.join(args.genbank_dir, "*.gb")) + glob.glob(os.path.join(args.genbank_dir, "*.gbk"))
    if not gb_paths:
        raise SystemExit("No GenBank files found.")

    gb_map = {os.path.basename(p): p for p in gb_paths}

    errors = []

    # if file column exists, use it; else best-effort: match by plasmid name if present
    file_col = None
    for c in df.columns:
        if c.lower() in ["filename", "genbank_file", "gb_file", "file"]:
            file_col = c
            break

    for i, row in df.iterrows():
        fp = None
        if file_col and isinstance(row[file_col], str) and row[file_col] in gb_map:
            fp = gb_map[row[file_col]]
        else:
            # heuristic: try IDs field pieces
            ids = str(row["IDs"]).replace(" ", "")
            # try match any filename containing first token
            token = ids.split(",")[0]
            candidates = [p for n,p in gb_map.items() if token in n]
            if candidates:
                fp = candidates[0]
        if not fp:
            errors.append(f"Row {i+1}: could not find GenBank file for IDs={row['IDs']}")
            continue

        rec = parse_genbank_relaxed(fp)
        seq = rec.seq

        # promoter check
        prom_ok = False
        prom_seq_ok = True
        for f in rec.features:
            if f.type == "promoter":
                lab = (f.qualifiers.get("label") or f.qualifiers.get("note") or [""])[0]
                lab = norm(lab)
                if row["Promoter_norm"] in lab:
                    prom_ok = True
                s = str(seq[int(f.location.start):int(f.location.end)]).upper()
                if len(s) != 35:
                    prom_seq_ok = False
        if not prom_ok:
            errors.append(f"{os.path.basename(fp)}: promoter label mismatch vs table ({row['Promoter_norm']})")
        if not prom_seq_ok:
            errors.append(f"{os.path.basename(fp)}: promoter sequence not 35 bp")

        # rbs check
        rbs_ok = False
        rbs_class_ok = False
        for f in rec.features:
            if f.type == "regulatory":
                rc = get_qual(f, "regulatory_class")
                if rc and rc.lower() == "ribosome_binding_site":
                    rbs_class_ok = True
                    lab = (f.qualifiers.get("label") or f.qualifiers.get("note") or [""])[0]
                    lab = norm(lab)
                    if row["RBS_norm"] in lab:
                        rbs_ok = True
        if not rbs_class_ok:
            errors.append(f"{os.path.basename(fp)}: missing regulatory_class=ribosome_binding_site for RBS")
        if not rbs_ok:
            errors.append(f"{os.path.basename(fp)}: RBS label mismatch vs table ({row['RBS_norm']})")

        # CDS internal stop check
        for f in rec.features:
            if f.type == "CDS":
                cds_seq = f.extract(seq)
                tr = str(cds_seq.translate(to_stop=False))
                if "*" in tr[:-1]:
                    errors.append(f"{os.path.basename(fp)}: CDS has internal stop codon(s)")

    if errors:
        print("VALIDATION FAILED:")
        for e in errors:
            print(" -", e)
        raise SystemExit(1)
    print("Validation OK: table and GenBank annotations are consistent.")
    return 0

if __name__ == "__main__":
    main()
