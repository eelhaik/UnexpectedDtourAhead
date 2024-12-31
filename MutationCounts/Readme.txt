# Pipeline for Generating Allele Counts Per Window  

This pipeline processes dbSNP (build 156) data to generate allele counts and visualize the percentage distribution of alleles across genomic windows.  

## Data Source  
The dbSNP data was downloaded from NCBI's FTP site (last accessed 9/9/2024):  
[https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/](https://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/)  
**File:** `GCF_000001405.40.vcf`  

## Workflow  

### 1. Extract Common Autosomal Variants (Excluding Indels)  
**Script:** `GetAutosomeAllelesNoIndelsCommon.txt`  
- Filters the dbSNP data to retain only common autosomal variants without indels.  
- **Output:** `output_no_indels_common.txt`  

### 2. Count Alleles Across Chromosomes  
**Script:** `CountBases.py`  
- Processes the filtered data to count alleles across chromosomes.  
- **Outputs:**  
  - `processed_output_all_chrs_50k.txt` (dbSNP allele counts)  
  - `processed_output_all_chrs_HGDP_50k.txt` (HGDP allele counts)  

### 3. Calculate Percentages Per Window  
- Input files:  
  - `processed_output_all_chrs_50k.txt`  
  - `processed_output_all_chrs_HGDP_50k.txt`  
- Generates allele counts for both dbSNP and HGDP datasets.  

### 4. Generate Plots  

#### Single-Chromosome Plot  
**Script:** `PlotWindowPercents.py`  
- Produces a percentage distribution plot for a single chromosome.  
- **Output:** `NC_000001.11_figure.png`  

#### Multi-Chromosome Plot  
**Script:** `PlotWindowPercents_all_chrs.py`  
- Produces percentage distribution plots for all chromosomes.  
- **Output:** `all_chromosomes_plot.png`  
