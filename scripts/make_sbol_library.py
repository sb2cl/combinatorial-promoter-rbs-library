#!/usr/bin/env python3
"""Generate SBOL3 combinatorial libraries (pSC101 and pGreen).

Inputs:
- metadata/library_table_from_paper.csv
- annotated GenBank files in genbank/

Outputs:
- sbol/library_pSC101.sbol3.xml
- sbol/library_pGreen.sbol3.xml

Usage:
  python scripts/make_sbol_library.py --base-uri "https://<your-namespace>/"
"""

import os, re, glob, argparse
from datetime import date

import pandas as pd
from Bio import SeqIO
from rdflib import Graph, Namespace, URIRef, Literal
from rdflib.namespace import RDF as RDFNS, DCTERMS, XSD

SBOL = Namespace("http://sbols.org/v3#")
SO   = Namespace("http://identifiers.org/so/")
SBO  = Namespace("http://identifiers.org/sbo/")
DNA_ENCODING = URIRef("http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html")

def get_qual(feat, key):
    vals = feat.qualifiers.get(key)
    return vals[0] if vals else None

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

def add_common(g, subject, display_id, name=None, desc=None):
    g.add((subject, SBOL.displayId, Literal(display_id)))
    if name:
        g.add((subject, SBOL.name, Literal(name)))
    if desc:
        g.add((subject, SBOL.description, Literal(desc)))

def make_sequence(g, uri, display_id, elements):
    g.add((uri, RDFNS.type, SBOL.Sequence))
    add_common(g, uri, display_id, name=display_id)
    g.add((uri, SBOL.elements, Literal(elements)))
    g.add((uri, SBOL.encoding, DNA_ENCODING))
    return uri

def make_component(g, uri, display_id, roles=None, seq_uri=None, desc=None):
    g.add((uri, RDFNS.type, SBOL.Component))
    add_common(g, uri, display_id, name=display_id, desc=desc)
    g.add((uri, SBOL.type, URIRef(SBO["0000251"])))  # DNA
    if roles:
        for r in roles:
            g.add((uri, SBOL.role, r))
    if seq_uri:
        g.add((uri, SBOL.sequence, seq_uri))
    return uri

def make_subcomponent(g, uri, display_id, instance_of, roles=None):
    g.add((uri, RDFNS.type, SBOL.SubComponent))
    add_common(g, uri, display_id, name=display_id)
    g.add((uri, SBOL.instanceOf, instance_of))
    if roles:
        for r in roles:
            g.add((uri, SBOL.role, r))
    return uri

def make_seq_constraint(g, uri, display_id, subj, obj):
    g.add((uri, RDFNS.type, SBOL.SequenceConstraint))
    add_common(g, uri, display_id)
    g.add((uri, SBOL.subject, subj))
    g.add((uri, SBOL.object, obj))
    g.add((uri, SBOL.restriction, SBOL.precedes))
    return uri

def make_variable_component(g, uri, display_id, variable_feature, variants):
    g.add((uri, RDFNS.type, SBOL.VariableComponent))
    add_common(g, uri, display_id)
    g.add((uri, SBOL.variable, variable_feature))
    for v in variants:
        g.add((uri, SBOL.variant, v))
    return uri

def make_comb_derivation(g, uri, display_id, template_comp, var_components):
    g.add((uri, RDFNS.type, SBOL.CombinatorialDerivation))
    add_common(g, uri, display_id, name=display_id, desc="Combinatorial promoter×RBS library encoded as SBOL3.")
    g.add((uri, SBOL.template, template_comp))
    g.add((uri, SBOL.strategy, SBOL.enumerate))
    for vc in var_components:
        g.add((uri, SBOL.variableComponent, vc))
    return uri

def extract_variant_sequences(genbank_paths):
    promoters = {}
    rbss = {}
    for fp in genbank_paths:
        rec = parse_genbank_relaxed(fp)
        seq = rec.seq
        for f in rec.features:
            if f.type == "promoter":
                lab = (f.qualifiers.get("label") or f.qualifiers.get("note") or ["promoter"])[0]
                m = re.search(r"(BBa[_-]J231\d+)", lab)
                if m:
                    lab = m.group(1).replace("-", "_")
                s = str(seq[int(f.location.start):int(f.location.end)]).upper()
                if len(s) == 35:
                    promoters.setdefault(lab, s)
            if f.type == "regulatory":
                rc = get_qual(f, "regulatory_class")
                if rc and rc.lower() == "ribosome_binding_site":
                    lab = (f.qualifiers.get("label") or f.qualifiers.get("note") or ["RBS"])[0]
                    m = re.search(r"(BBa[_-](?:B\d{4}|J\d{5}))", lab)
                    if m:
                        lab = m.group(1).replace("-", "_")
                    s = str(seq[int(f.location.start):int(f.location.end)]).upper()
                    rbss.setdefault(lab, s)
    return promoters, rbss

def build_library_sbol(backbone_name, promoter_set, rbs_set, base_uri):
    g = Graph()
    g.bind("sbol", SBOL); g.bind("dcterms", DCTERMS); g.bind("so", SO); g.bind("sbo", SBO)

    doc_uri = URIRef(f"{base_uri}{backbone_name}/")
    g.add((doc_uri, RDFNS.type, SBOL.Document))
    g.add((doc_uri, DCTERMS.created, Literal(str(date.today()), datatype=XSD.date)))
    g.add((doc_uri, DCTERMS.title, Literal(f"pIAKA combinatorial library ({backbone_name}) - SBOL3")))

    promoter_comps = []
    for lab, seq in sorted(promoter_set.items()):
        seq_uri = URIRef(f"{base_uri}seq/{lab}")
        comp_uri = URIRef(f"{base_uri}component/{lab}")
        make_sequence(g, seq_uri, f"seq_{lab}", seq)
        make_component(g, comp_uri, lab, roles=[URIRef(SO["0000167"])], seq_uri=seq_uri,
                       desc="Constitutive promoter (BioBrick J231xx family).")
        promoter_comps.append(comp_uri)

    rbs_comps = []
    for lab, seq in sorted(rbs_set.items()):
        seq_uri = URIRef(f"{base_uri}seq/{lab}")
        comp_uri = URIRef(f"{base_uri}component/{lab}")
        make_sequence(g, seq_uri, f"seq_{lab}", seq)
        make_component(g, comp_uri, lab, roles=[URIRef(SO["0000139"])], seq_uri=seq_uri,
                       desc="Ribosome binding site (RBS).")
        rbs_comps.append(comp_uri)

    gfp_uri = URIRef(f"{base_uri}component/GFPmut3")
    make_component(g, gfp_uri, "GFPmut3", roles=[URIRef(SO["0000316"])], desc="Reporter coding sequence (GFPmut3).")
    term_uri = URIRef(f"{base_uri}component/Terminator")
    make_component(g, term_uri, "Terminator", roles=[URIRef(SO["0000141"])], desc="Transcriptional terminator.")
    bb_uri = URIRef(f"{base_uri}component/{backbone_name}_backbone")
    make_component(g, bb_uri, f"{backbone_name}_backbone", roles=[URIRef(SO["0000755"])], desc=f"Plasmid backbone ({backbone_name}).")

    template_uri = URIRef(f"{base_uri}component/ExpressionCassette_{backbone_name}")
    make_component(g, template_uri, f"ExpressionCassette_{backbone_name}", desc=f"Template construct in {backbone_name} backbone.")

    sc_prom = URIRef(f"{base_uri}feature/{backbone_name}/promoter")
    sc_rbs  = URIRef(f"{base_uri}feature/{backbone_name}/rbs")
    sc_gfp  = URIRef(f"{base_uri}feature/{backbone_name}/gfp")
    sc_term = URIRef(f"{base_uri}feature/{backbone_name}/terminator")
    sc_bb   = URIRef(f"{base_uri}feature/{backbone_name}/backbone")

    make_subcomponent(g, sc_prom, f"{backbone_name}_promoter", promoter_comps[0], roles=[URIRef(SO["0000167"])])
    make_subcomponent(g, sc_rbs,  f"{backbone_name}_rbs",      rbs_comps[0], roles=[URIRef(SO["0000139"])])
    make_subcomponent(g, sc_gfp,  f"{backbone_name}_gfp",      gfp_uri, roles=[URIRef(SO["0000316"])])
    make_subcomponent(g, sc_term, f"{backbone_name}_terminator", term_uri, roles=[URIRef(SO["0000141"])])
    make_subcomponent(g, sc_bb,   f"{backbone_name}_backbone_sc", bb_uri, roles=[URIRef(SO["0000755"])])

    for feat in [sc_prom, sc_rbs, sc_gfp, sc_term, sc_bb]:
        g.add((template_uri, SBOL.feature, feat))

    make_seq_constraint(g, URIRef(f"{base_uri}constraint/{backbone_name}/p_before_r"), f"{backbone_name}_p_before_r", sc_prom, sc_rbs)
    make_seq_constraint(g, URIRef(f"{base_uri}constraint/{backbone_name}/r_before_g"), f"{backbone_name}_r_before_g", sc_rbs, sc_gfp)
    make_seq_constraint(g, URIRef(f"{base_uri}constraint/{backbone_name}/g_before_t"), f"{backbone_name}_g_before_t", sc_gfp, sc_term)

    vc_prom = URIRef(f"{base_uri}variable/{backbone_name}/promoter")
    vc_rbs  = URIRef(f"{base_uri}variable/{backbone_name}/rbs")
    make_variable_component(g, vc_prom, f"{backbone_name}_var_promoter", sc_prom, promoter_comps)
    make_variable_component(g, vc_rbs,  f"{backbone_name}_var_rbs",      sc_rbs,  rbs_comps)

    cd_uri = URIRef(f"{base_uri}derivation/Library_{backbone_name}")
    make_comb_derivation(g, cd_uri, f"Library_{backbone_name}", template_uri, [vc_prom, vc_rbs])
    return g

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-uri", default="https://example.org/pIAKA/", help="SBOL namespace base URI")
    ap.add_argument("--table", default="metadata/library_table_from_paper.csv", help="CSV from paper table")
    ap.add_argument("--genbank-dir", default="genbank", help="Directory containing .gb files")
    ap.add_argument("--outdir", default="sbol", help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.table)

    gb_paths = glob.glob(os.path.join(args.genbank_dir, "*.gb")) + glob.glob(os.path.join(args.genbank_dir, "*.gbk"))
    if not gb_paths:
        raise SystemExit(f"No GenBank files found in {args.genbank_dir}/")

    promoters, rbss = extract_variant_sequences(gb_paths)

    def norm(x): 
        return str(x).strip().replace("-", "_")

    df["Promoter_norm"] = df["Promoter"].map(norm)
    df["RBS_norm"] = df["RBS"].map(norm)

    df_psc = df[df["ORI"].str.contains("pSC101", case=False, na=False)]
    df_pgr = df[df["ORI"].str.contains("pGreen", case=False, na=False)]

    prom_psc = {k: promoters[k] for k in sorted(set(df_psc["Promoter_norm"])) if k in promoters}
    rbs_psc  = {k: rbss[k] for k in sorted(set(df_psc["RBS_norm"])) if k in rbss}

    prom_pgr = {k: promoters[k] for k in sorted(set(df_pgr["Promoter_norm"])) if k in promoters}
    rbs_pgr  = {k: rbss[k] for k in sorted(set(df_pgr["RBS_norm"])) if k in rbss}

    g1 = build_library_sbol("pSC101", prom_psc, rbs_psc, args.base_uri)
    g2 = build_library_sbol("pGreen", prom_pgr, rbs_pgr, args.base_uri)

    out1 = os.path.join(args.outdir, "library_pSC101.sbol3.xml")
    out2 = os.path.join(args.outdir, "library_pGreen.sbol3.xml")
    g1.serialize(destination=out1, format="xml")
    g2.serialize(destination=out2, format="xml")

    print("Wrote:", out1)
    print("Wrote:", out2)

if __name__ == "__main__":
    main()
