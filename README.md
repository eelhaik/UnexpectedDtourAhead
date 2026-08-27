# Version 2 Scripts

This folder contains two independent pieces of code used in the paper:

1. **`Example C++ code.txt`** — parses the genomic input files (`MASTER_CHR*` files) and computes windowed heterozygosity, ABBA/BABA counts, D statistics, and a relative mutation-rate estimator.
2. **`ESM10.zip`** — a SLiM forward simulation package (population-genetic simulation), together with the triplet mutation-rate matrix `mm_forSLiM.csv` that it depends on.
3. **1KG VCF Files** -Compressed files of the type MASTER_CHR*.txt are available here: https://zenodo.org/records/22119343

These are examples/starting points, not turnkey pipelines — both read-me files below stress that the user is expected to adapt the code to the specific output they want.

---

## 1. C++ analysis code

**Files:** `Example C++ code.txt`, `C++ example code read me.txt`

### What it does

Reads per-chromosome master files (biallelic SNPs from the high-coverage 1000 Genomes data, scored against Vindija Neanderthal, Altai Neanderthal, Denisovan, and chimpanzee outgroups), and for user-defined genomic windows:

* Computes the heterozygosity difference between Africans and non-Africans, bins it onto a 0–99 scale.
* Computes ABBA/BABA counts and D statistics for 3 site classes (XBXA, BBBA, BBAA) × 3 population comparisons (EUR–AFR, EAS–AFR, EAS–EUR).
* Computes a relative mutation-rate estimator per heterozygosity bin (Europe vs. Africa).
* Writes everything to a single tab-delimited output file, one row per heterozygosity bin.

### Input file format

Each `MASTER_CHR<N>.txt` file (one per chromosome) is whitespace-delimited with a header row, and per-site columns in this order:

```text
CHR  LOC  REF  ALT  fqGBR..fqPEL (26 population alt-allele frequencies)  fqAFR  fqASIA  MAJ1  MAJ2  TRIP  OGP
```

* `fqAFR`, `fqASIA` — balanced-sample alt-allele frequencies for Africa / outside Africa.
* `TRIP` — mutating ancestral triplet code (1–128, see field description in `C++ example code read me.txt`).
* `OGP`/`CP` — outgroup conformation code, with taxa ordered as **Vindija Neanderthal, Altai Neanderthal, Denisovan, chimpanzee**: 0=X1X0, 1=X0X1, 2=0000, 3=1111, 4=0001, 5=1110, 6=1100, 7=0011.

The code expects files named `MASTER_CHR<chromosome>.txt` in a single input directory.

### How to run

1. Open a C++ compiler/IDE (the authors used **Visual Studio 2022**; any C++17-capable compiler with `<filesystem>` support will work).
2. Create a new console project and paste the contents of `Example C++ code.txt` into `main.cpp`.
3. Build/compile (C++17 or later).
4. Run the executable. It will prompt interactively for:

   * **Input directory** containing the `MASTER_CHR*.txt` files.
   * **Output filename** (include the full path, or it is written to the current directory).
   * **First chromosome** and **last chromosome** (default 1–22).
   * **Genomic window size** (default 50,000 bp).
5. The program processes each chromosome in turn and writes a single results file containing, per heterozygosity-difference bin: raw ABBA/BABA counts, D statistics, and the mutation-rate estimator (`MU_AFR`, `MU_NAF`, `RelMU`).

### Notes

* This is explicitly an **example/starting point** — it correctly parses the input and reproduces one specific analysis (heterozygosity-partitioned D and relative mutation rate). Other outputs described in the paper require editing the code (the authors suggest this is straightforward with AI assistance, or on request).
* If you change the window size, heterozygosity differences scale proportionally — you may need to adjust the binning/scaling logic (see `heterozygosity_bin()`) to keep values in a sensible range.
* Genome-wide analysis is recommended when partitioning by triplet, since per-triplet ABBA/BABA counts are small and noisy.

---

## 2. SLiM simulation package (`ESM10.zip`)

**Files inside the zip:** `NUC_GENERAL.txt` (SLiM simulation script), `SLiM_read.me.txt` (instructions).

**Companion file (already in this folder):** `mm_forSLiM.csv` — the 64×4 triplet-specific mutation-rate matrix used when the full triplet mutation model is selected.

### What it does

`NUC_GENERAL.txt` is a nucleotide-based SLiM model simulating a chimpanzee outgroup and a Hominin lineage that splits into Neanderthal and African populations, with a non-African population subsequently splitting off and receiving Neanderthal introgression. It simulates 1 chimpanzee, 50 Africans, 1 Neanderthal, and 50 non-Africans (all diploid) over a 500 kb sequence, and outputs per-site position, reference/alternate bases, mutation provenance, and per-population allele frequencies, plus final heterozygosity values for Africans and non-Africans.

### Requirements

* **SLiM** version 5.1 (or later) — https://messerlab.org/slim/
* Designed to run under **WSL (Windows Subsystem for Linux)** using a bash wrapper script, but any Linux/macOS shell with SLiM installed will work.
* 500 kb ancestral sequence files in FASTA format, named `segmentAA_<N>.fa`, sampled from human chromosome 1 and located in a path configured by the user (see below).
* If using the full triplet mutation model (`FLAG=2`), `mm_forSLiM.csv` must be present in the working directory from which SLiM is launched.

We used **500 different initializing segments from human chromosome 1**. For seed numbers 501–1000, these sequences are recycled (501 → 1, 502 → 2, etc.). Users requiring more replicates, or who prefer a unique sequence in every run, can easily extend the initialization library.

### Setup

1. Extract `ESM10.zip` (contains `NUC_GENERAL.txt` and `SLiM_read.me.txt`).
2. Place `mm_forSLiM.csv` alongside `NUC_GENERAL.txt` (needed only for `FLAG=2`).
3. Edit the `PATH` constant near the top of `NUC_GENERAL.txt` (currently `/mnt/h/SimsN/`) to point to the directory containing your `segmentAA_<N>.fa` sequence files.

### How to run

Run from a bash shell (WSL or Linux) with SLiM on the PATH. Example — runs 100 replicates:

```bash
a=100
for k in $(seq $((a-99)) $a); do
    slim -d SEED=$k -d NEA=5 -d INCR=1 -d MIG=0.005 -d DROP=0.0 -d FLAG=2 NUC_GENERAL.txt
done
```

**Parameters** (passed via `-d`):

| Flag   | Meaning                                                                                                                                                                 |
| ------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SEED` | Replicate number; selects the corresponding `segmentAA_<SEED>.fa` ancestral sequence. For seeds >500, the 500 initialization sequences are recycled as described above. |
| `NEA`  | Neanderthal population size, in thousands (paper default: `5`).                                                                                                         |
| `INCR` | Rate of rise of the general human mutation rate (`1` = no rise; max sensible value `5`).                                                                                |
| `MIG`  | Introgression (migration) rate from Neanderthals into non-Africans (~`0.005`–`0.01` gives realistic D ≈ 0.05).                                                          |
| `DROP` | Proportional reduction in non-African mutation rate (e.g. `0.25` = 25% reduction).                                                                                      |
| `FLAG` | Mutation model: `1` = simple transition/transversion (Kimura) model; `2` = full triplet model (requires `mm_forSLiM.csv`).                                              |

Each run produces an output file (one per `SEED`) with per-site position, reference base, alternate base(s) tagged with generation and population of origin, per-population allele frequencies (p1=chimpanzee, p2=Africa, p3=Neanderthal, p4=non-Africa), and Africa/non-Africa heterozygosity at the end.

### Notes

* The authors also used adapted variants of the script with a non-African population size of 10,000 and either no bottleneck or a step bottleneck. These variants are not included here but can be easily derived by editing the demography section of `NUC_GENERAL.txt`, **or can be requested from the authors**.
* `mm_forSLiM.csv` has 64 rows (one per ancestral triplet, central-base ordering as used in the C++ code) × 4 columns (per-target-base relative mutation rate).
