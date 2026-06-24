# EF-Tu Structure-Guided Clustering and Phylogenetic/Ecological Analysis

Code accompanying:

**D81 Mutants Reveal Hidden EF-Tu Diversity while Natural Sequences Preserve Aspartate**
Jordan L. Johnson and Yuhong Wang, *Computational and Structural Biotechnology Journal* (accepted, 2026).

This repository implements a structure-guided sequence-analysis pipeline that uses AlphaFold2/ColabFold predictions to (1) generate diverse predicted conformations from wild-type and mutant query sequences, (2) sort those structures with a two-step RMSD clustering strategy, and (3) link each conformational cluster back to its contributing sequences and their protein identity, taxonomy, and ecology. The same framework is applied to EF-Tu D81 mutants (D81A, D81F, D81K) versus wild type, but is general to any conserved protein.

---

## 1. Overview of the pipeline

The workflow proceeds in three stages, mirroring Figure 2 of the paper. The numbers in parentheses below correspond to the cell-number comments inside the notebook (`#1`, `#2`, …), which are also the numbers cited in Figure 2.

**Stage 1 — Generate high-confidence structures (cells #1–#5)**
A query sequence is folded with ColabFold to produce an MSA (`.a3m`); the MSA is clustered (DBSCAN, adapted from AF-Cluster); clusters above a size threshold are folded; and only predictions with mean pLDDT > 70 are retained.

**Stage 2 — Two-step conformational clustering (cells #6–#10)**
All retained structures from all four queries are pooled and source-tagged. Structures are first globally aligned over the rigid domains II–III (residues 212–394), then clustered by pairwise RMSD of the domain I terminal helix (residues 183–199). A second, finer sorting is applied over the switch I loop (residues 54–58). Clusters are read off a dendrogram / off-diagonal RMSD trace.

**Stage 3 — Sequence, taxonomic, and ecological annotation (cells #11–#28)**
Sequences contributing to the retained structures are pooled, deduplicated by UniRef100 ID (with all contributing query sources consolidated into the label), annotated with protein identity and taxonomy fetched from UniProt, used to build an IQ-TREE phylogeny with ancestral reconstruction, and visualized with iTOL. Downstream cells generate the LOGO, species, genus, and set-overlap analyses.

---

## 2. Requirements

### Python environment (conda)

Two environment files are provided:

- `environment_full.yml` — the complete, exact environment (Linux, with build strings) used to produce the published results. Use this for byte-for-byte reproduction on Linux.
- `environment.yml` — the same packages without build strings, for cross-platform installation.

Create and activate the environment:

```bash
conda env create -f environment.yml
conda activate myenv
```

Key Python packages: `python=3.10`, `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `biopython`, `ete3`, `pymol-open-source`, `fastcluster`, `logomaker`, `matplotlib-venn`, `upsetplot`, `requests`, `openpyxl`.

### External tools (installed separately, not in the conda file)

These are called by the notebook but are not Python packages, so install them yourself and update the paths in cell `#1`:

| Tool | Used in | Notes |
|------|---------|-------|
| **ColabFold (localcolabfold)** v1.5.x | cells #2, #4 | provides `colabfold_batch`; set its path in `COLAB` (cell #1) |
| **AF-Cluster scripts** (`ClusterMSA_min3.py`, `utils.py`) | cell #3 | adapted from Wayment-Steele et al. 2024; place in `default_path` |
| **AlphaFold3 server** | (separate, for Figure 1/S7–S13) | run at https://alphafoldserver.com/ ; inputs listed in the paper Methods |
| **IQ-TREE 2** | cells #20, #21 | must be on `PATH` as `iqtree2` |
| **iTOL** (web) | visualization | tree + datasets uploaded to https://itol.embl.de/ |

---

## 3. Repository structure

```
.
├── README.md
├── myenv_Tu.yml         		  # exact Linux environment
├── Tu_2025_code.ipynb           # the full pipeline notebook
├── ClusterMSA_min3.py           # MSA clustering (from AF-Cluster, adapted)
├── utils.py                     # helper functions for clustering
└── (working directories created at runtime: SEQ/, Cluster/, Combined_sel_pdb/, tree/ …)
```

At runtime the notebook creates a working tree under `default_path`:

```
Tu_paper2025/
├── SEQ/                         # input query .a3m (one at a time)
├── Cluster/                     # clustered a3m
│   └── folding/                 # folded pdbs + pLDDT.csv
│       └── sel_pdb/             # pLDDT>70 structures, source-tagged
└── Combined_sel_pdb/            # pooled structures for clustering
    └── 1st_212-394_RMSD/        # global alignment + RMSD
        └── 2nd_183-199_RMSD/    # helix clustering
            └── tree/            # combined summary, fasta, IQ-TREE, iTOL datasets
```

---

## 4. How to run

The notebook is meant to be run **cell by cell, in order**, because each stage produces files consumed by the next. A few stages (folding, IQ-TREE) are long-running and a couple require a kernel restart (noted in the cell comments).

1. **Set paths and parameters (cell #1).** This is the single control cell. Edit:
   - `default_path` — your project root
   - `COLAB` — path to `colabfold_batch`
   - thresholds: cluster size (`size_threshold`, default 4 KB ≈ ≥7 sequences), `pLDDT_threshold` (default 70), the alignment range (`212-394`) and the clustering range (`183-199`).

2. **Fold and cluster (cells #2–#5).** Place one query `.a3m` in `SEQ/`, fold, cluster the MSA, fold the clusters, and keep pLDDT > 70 structures. Repeat for each query (WT, D81A, D81F, D81K), changing the source tag.

3. **Pool and two-step cluster (cells #6–#10).** Combine all source-tagged structures, align over domains II–III, compute pairwise RMSD over helix 183–199, and read clusters from the dendrogram / off-diagonal trace. Cell #7 may need a kernel restart and a rerun of cell #1 first (noted in the comment).

4. Annotate and analyze (cells #11–#28).
Pool and align by UniRef ID (cells #11–#13). The selected .a3m files from all four queries are gathered and combined cluster by cluster (cell #11), then written into a columnized summary (cell #12). Sequences are keyed by their UniRef100 ID. When the same UniRef100 ID appears in more than one query's MSA, it is collapsed into a single representative sequence, and every contributing source is accumulated into that entry's tag rather than duplicated as separate rows (cell #13). This is what produces labels like >1_WT_D81A_D81K_<UniRefID>, where the leading number is the conformational cluster and the middle tokens record exactly which queries recruited that sequence.
Retrieve the aligned sequences (cell #14). Using the deduplicated ID list, the corresponding aligned sequences are pulled back out of the combined FASTA so that downstream analyses operate on one sequence per unique ID.
Fetch protein identity and taxonomy from UniProt (cells #15–#19). Only after deduplication, each unique UniRef100 ID is queried against UniProt to retrieve its protein name and taxonomic lineage (cell #15–#16), producing the three-column protein / taxon / UniRef100 table. A short-name mapping is applied so protein identities can be grouped consistently (cells #17–#19). The final unified label therefore encodes cluster, query source(s), UniRef100 ID, protein, and organism in a single header.
Build the phylogeny (cells #20–#21). The deduplicated, labeled alignment is used to construct a maximum-likelihood tree with IQ-TREE 2 (ModelFinder, 1000 ultrafast bootstraps), followed by ancestral sequence reconstruction. Cell #21 can resume the ancestral step alone if the combined run is interrupted.
Visualize and analyze (cells #22–#28). iTOL annotation datasets are generated for protein identity, conformational cluster, and species rings (cells #22–#24); position-81 LOGO plots are made for inclusive/exclusive query subsets (cell #25); and species- and genus-level distributions, ID extraction, and set overlaps are computed (cells #26–#28).

---

## 5. Reproducing the paper figures

| Figure | Produced by |
|--------|-------------|
| Fig. 2 (workflow) | schematic; cell numbers match the `#N` comments |
| Fig. 3A (dendrogram / RMSD heatmap + inset) | cell #9 (RMSD) and cell #10 (dendrogram + off-diagonal trace) |
| Fig. 3B (representative structures) | PyMOL sessions from cells #7–#8 |
| Fig. 4A (phylogenetic tree + rings) | IQ-TREE (cells #20–21) + iTOL datasets (cells #22–24) |
| Fig. 4B–C (species analyses) | cells #23, #26, #27 |
| Fig. 5 (LOGO of position 81) | cell #25 |
| Fig. S11 (switch I sorting) | second-range sorting, cell #8 |

---

## 6. Key parameters (all set in cell #1)

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `pLDDT_threshold` | 70 | minimum mean pLDDT to keep a predicted structure |
| `size_threshold` | 4 KB | minimum MSA cluster size to fold (≈ ≥7 sequences) |
| `residue_range1` | `212-394` | domains II–III, used as the alignment reference |
| (helix range) | `183-199` | domain I terminal helix, used for global clustering |
| (switch I range) | `54-58` | switch I loop, used for second-step local clustering |
| RMSD cluster cutoff | 10 Å | threshold on the off-diagonal RMSD trace |

---

## 7. Data availability

Query sequences are E. coli EF-Tu (UniProt P0CE48) and its D81A/D81F/D81K variants. Retrieved MSAs, predicted structures, and the final annotated summary tables are generated by the pipeline; intermediate large files are not stored in the repository but are fully reproducible from the steps above.

---

## 8. Notes and caveats

- **AlphaFold predictions are MSA-dependent**, not free-energy samples; conformational diversity here reflects MSA subsampling, and confidence metrics (pLDDT, PAE) report internal consistency, not physical accuracy.
- **MSA composition is query-dependent** (query sequence, database version, search depth); the recruited sequence set is not an unbiased census of natural diversity.
- **Some steps are not fully scriptable**: AlphaFold3 modeling is run on the public server, and final tree styling is done in iTOL.
- **Cell numbers** in this README and in Figure 2 refer to the `#N` comments in the notebook, not to Jupyter execution counts (which reset on rerun).

---

## 9. Citation

If you use this code, please cite the paper:

> Johnson, J.L., and Wang, Y. (2026). D81 Mutants Reveal Hidden EF-Tu Diversity while Natural Sequences Preserve Aspartate. *Computational and Structural Biotechnology Journal*.

and the key tools it builds on: ColabFold (Mirdita et al. 2022), AF-Cluster (Wayment-Steele et al. 2024), AlphaFold2 (Jumper et al. 2021) and AlphaFold3 (Abramson et al. 2024), IQ-TREE 2 (Minh et al. 2020), and iTOL v6 (Letunic & Bork 2024).
