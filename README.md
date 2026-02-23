## SUMOylation

This repository contains scripts for the data analysis of the results from the study for developing MaxSBM. The work presents MaxSBM, a module in MaxQuant, developed for improved site-specific identification of
SUMOylation sites, using mass-spectrometry-based proteomics.

Documentation of the module can be found here: https://cox-labs.github.io/coxdocs/MaxSBM.html

### The work is presented in;

- **Lennartsson C, Kyriakidou P, Nielsen ML, Olsen JV, Cox J, Hendriks IA**  
  *Improved peptide search for identification of SUMO and sequence-based modifications, in MaxSBM.*  
  bioRxiv (2025).  
  https://www.biorxiv.org/content/10.1101/2025.08.27.672604v1.abstract

### Datasets

#### SUMO-HEK (human cell lines)
- **Hendriks IA, Lyon D, Su D, Skotte NH, Daniel JA, Jensen LJ, Nielsen ML** (2018)  
  *Site-specific characterization of endogenous SUMOylation across species and organs.*  
  *Nature Communications*, **9**(1), 2456.  
  https://doi.org/10.1038/s41467-018-04957-4

#### SUMO-MEC (mouse embryonic cells)
- **Theurillat I, Hendriks IA, Cossec J-C, Andrieux A, Nielsen ML, Dejean A** (2020)  
  *Extensive SUMO modification of repressive chromatin factors distinguishes pluripotent from somatic cells.*  
  *Cell Reports*, **32**(11), 108146.  
  https://doi.org/10.1016/j.celrep.2020.108146

#### SUMO-Adip (mouse adipocytes)
- **Zhao X, Hendriks IA, Le Gras S, Ye T, Ramos-Alonso L, Nguéa PA, et al.** (2022)  
  *Waves of SUMOylation support transcription dynamics during adipocyte differentiation.*  
  *Nucleic Acids Research*, **50**(3), 1351–1369.  
  https://doi.org/10.1093/nar/gkac027

#### Ub-LysC (human cell lines)
- **Akimov V, Barrio-Hernandez I, Hansen SVF, Hallenborg P, Pedersen AK, et al.** (2018)  
  *UbiSite approach for comprehensive mapping of lysine and N-terminal ubiquitination sites.*  
  *Nature Structural & Molecular Biology*, **25**, 631–640.  
  https://doi.org/10.1038/s41594-018-0084-y

### Notebooks and scripts

This repository contains the following notebooks and scripts:

#### Method development and ion characterization
- **[Diagnostic ion mining](Diagnostic_peak_mining.ipynb)**  
  Identification and characterization of diagnostic fragment ions associated with SUMO.

- **[d-ion series in search](d_ion_series_in_search.ipynb)**  
  Assessment of the impact of including d-ion series during peptide database searching.

- **[p-ion optimisations (SUMO)](p_ion_optimisations_sumo.ipynb)**  
  Assessment of the impact of including p-ion series during peptide database searching.

- **[Which p-SUMO ions](which_p_sumo.ipynb)**  
  Evaluation of p-ion inlcusions and their contribution to peptide identification.

#### Spectral inspection and visualization
- **[Spectral viewer](spectral_viewer.ipynb)**  
  Interactive visualization of MS/MS spectra for inspection of fragment ion series and modification-specific peaks.

- **[Benchmark plots](theoretical_peptide.ipynb)**  
  Generation and analysis of theoretical peptide and fragment ion compositions for SUMoylated peptides. Used for problem formulation.

#### Dataset-specific analyses of SUMOylatated peptides
- **[Main SUMO-HEK analysis](Main_result_SUMO_HEK.ipynb)**  
  Primary analysis and summary of SUMOylation results in human HEK cell lines.

- **[Mouse embryonic cells analysis](mouce_embryonic_cells_analysis.ipynb)**  
  Analysis of SUMOylation data from mouse embryonic cells (MEC).

- **[Mouse adipocytes analysis](mouse_adipocytes_analysis.ipynb)**  
  Analysis of SUMOylation data from mouse adipocytes.

#### Protein-level analysis
- **[Protein groups (SUMO)](protein_groups_sumo.ipynb)**  
  Protein-level summarization and grouping of SUMO-modified peptides.

- **[Protein groups (adipocytes)](protein_groups_adiposites.ipynb)**  
  Protein group analysis for adipocyte datasets.

- **[Protein groups (MEC)](protein_groups_mec.ipynb)**  
  Protein group analysis for mouse embryonic fibroblast datasets.

- **[Protein inference](protein_inference.py)**  
  Python script implementing protein inference from uniprot metadata.

#### Data processing and utilities
- **[Search output processing](Process_dataframes.ipynb)**  
  Post-processing, filtering, and restructuring of search engine output tables.

- **[Percolator converter](percolator_converter.ipynb)**  
  Conversion of MaxQuant search results into Percolator-compatible input formats.
