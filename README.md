[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.ZENODO_RECORD_ID.svg)](https://doi.org/10.5281/zenodo.ZENODO_RECORD_ID)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![validate-library](https://github.com/sb2cl/combinatorial-promoter-rbs-library/actions/workflows/validate-library.yml/badge.svg)](https://github.com/sb2cl/combinatorial-promoter-rbs-library/actions/workflows/validate-library.yml)
[![SBOL3](https://img.shields.io/badge/SBOL-3-blue.svg)](sbol/)

# Combinatorial promoter–RBS library

This repository contains:
- Annotated GenBank files (`genbank/`)
- SBOL3 design files with *formal combinatorial derivations* split by backbone (`sbol/`)
  - `library_pSC101.sbol3.xml`
  - `library_pGreen.sbol3.xml`

## SBOL3 base URI
The SBOL3 files were generated with a placeholder base namespace:

    https://github.com/sb2cl/combinatorial-promoter-rbs-library/

Before publishing, you may wish to replace this with your final repository namespace
(e.g., your GitHub Pages URL or another stable URI). A simple global search/replace
is sufficient.

## Suggested citation
When you create a Zenodo release, cite the Zenodo DOI here and in the manuscript.

## Reproducibility

### Validate the library
```bash
python scripts/validate_library.py --table metadata/library_table_from_paper.csv --genbank-dir genbank
```

### Regenerate SBOL3 files
```bash
python scripts/make_sbol_library.py --base-uri "https://<your-namespace>/" --table metadata/library_table_from_paper.csv --genbank-dir genbank --outdir sbol
```

The GitHub Action `.github/workflows/validate-library.yml` runs these checks automatically on every push/PR.

