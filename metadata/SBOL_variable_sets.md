## SBOL3 combinatorial design (recommended)

Define a template Component (expression cassette) with ordered SubComponents:
Promoter -> RBS -> GFPmut3 -> Terminator.

Then define a CombinatorialDerivation with:
- VariableComponent `vc_promoter` targeting the Promoter SubComponent:
  variants = {BBa_J23106, BBa_J23102, BBa_J23101, BBa_J23100}
- VariableComponent `vc_rbs` targeting the RBS SubComponent:
  variants = {BBa_B0030, BBa_B0032, BBa_J61100, BBa_J61101, BBa_B0034}

Optionally, model backbone choice as:
- either separate top-level Components for pSC101 and pGreen backbones (recommended for clarity),
- or a third VariableComponent `vc_backbone` if you want a single derivation spanning both.

The file `metadata/library_table_from_paper.csv` provides the authoritative mapping of construct IDs to
Promoter/RBS/backbone combinations.
