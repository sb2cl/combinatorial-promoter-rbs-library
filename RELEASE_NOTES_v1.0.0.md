# Release v1.0.0

## Combinatorial Promoter–RBS Expression Library

This release provides the fully annotated DNA sequences and formal SBOL3 specification of a combinatorial promoter–RBS expression library constructed in pSC101 and pGreen backbones.

### Contents

- Annotated GenBank files for all plasmid constructs
- SBOL3 combinatorial derivations encoding the promoter × RBS design space
- Machine-readable metadata linking construct identifiers to biological parts
- Reproducibility scripts and automated validation workflow

### Design Space

The library encodes combinations of:

- Constitutive sigma70 promoters (BioBrick J231xx family)
- Ribosome binding sites (B003x and J611xx series)
- GFPmut3 reporter
- Standard transcriptional terminator
- Two backbone contexts: pSC101 and pGreen

### Reproducibility

The repository includes:

- `scripts/make_sbol_library.py` to regenerate SBOL3 derivations
- `scripts/validate_library.py` to validate GenBank annotations against metadata
- Continuous integration checks via GitHub Actions

### Citation

Please cite the associated publication and this repository DOI (Zenodo).

---

SB2CL — Synthetic Biology and Biosystems Control Lab
Universitat Politècnica de València
