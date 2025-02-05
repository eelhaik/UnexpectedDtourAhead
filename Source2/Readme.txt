**Instructions for Using Source2.cpp**

**1. Introduction**

This program calculates aspects of the introgression statistic D, broken down into transitions and transversions (TS/TV), which 3-base motif (triplet) has mutated and the frequencies of the alleles in different populations. 
Data are analysed both in non-overlapping genomic windows of fixed size and also as blocks of 100,000 contiguous qualifying sites, the former to look at variation across the genome, the latter to generate confidence intervals
The output generated is a tab delimited file that is best viewed by opening in Excel or equivalent and starts with a summary of the parameter values supplied by the user.
Below the parameter values are a series of Tables corresponding to specific Figures in the paper. Fuller explanations for each are in the text of the paper.
Note, analysing just the smallest chromosome 22 results in some outputs containing low or zero counts for the rarest base combinations.
Note also, most analyse are conducted with window size of 50Kb. Where this is not the case, the program will need to be run again using a different size.

The first block contains data for Figure 4. C_dist is the branch length asymmetry to the chimpanzee. N>P1 and N>P2 are the binned branch lengths to the Altai Neanderthal. D'N is the Neanderthal branch length asymmetry. For each bin we also calculate ABBA, BABA and D.

The second block contains data for Figure 5. Here, three taxon combinations are explored: comparisons within P1; comparisons within P1 where the second individual is spiked with Neanderthal alleles; comparisons between P1 and P2.
The output gives raw counts for each branch length (first 6 columns), followed by the corresponding D' values, in order. 

The third block contains the same data as the second block, but here the data are binned on a much finer scale around zero.

The forth block contains data for Figure 6. Here, D is calculated for the three different taxon combinations in block 2, with separate values calculated for different classes of site: transitions and transversions and mutating triplet. 

The fifth block contains data for the righthand panel in Figure 11. Here, D is calculated in bins based on the slope of the relationship between heterozygosity across Eurasia and distance from Africa. D is calculated both overall, and separately depending on whether the human major allele is the chimpanzee A or the Neanderthal B.

The sixth block contains data for the lefthand panel in Figure 11. Here, the three versions of D are calculated as in the previous block, but here the binning is based on heterozygosity difference between Bantu and Japanese.

The seventh and last block contains data for Figure S3. Here, data are binned by the human to chimpanzee branch asymmetry and used to calculate the human Neanderthal branch asymmetry.
Three different Neanderthal branch length asymmetries are calculated, depending on the based carried by the non-humans: chimpanzee=A, Neanderthal=A; chimpanzee=A, Neanderthal = B; chimpanzee=B, Neanderthal=A. 


**2. Prerequisites**

*   The program is designed to be compiled and run in a C++ environment.
*   The following input files should be located in the same directory:
    *   HGDP VCF files (e.g., `hgdp.v0.5.archaics.chr22.vcf`)
    *   Reference FASTA files (e.g., `chr22.fa`)
    *   Sample information file (`inpopsHGDP.txt`)
    *   Population characterists file ('HGDP_pop_characters.txt')

*   Filenames must  be exact and are case-sensitive.

**3. Initial Setup**

*   The program uses a **graphical interface** to select the directory with the input files. The program will prompt you to select a folder using a dialog box.
*   After selecting the folder, the program will load the population data from `inpopsHGDP.txt` and the reference sequence data from the `chr[chromosome number].fa` files.
*   The program will then prompt you for an **output file name**. The output file will be a `.txt` file and is saved into the selected directory. If '.txt' is not in the filename it will be added automatically.

* In Linux, there is no graphic interface. To compile the code do:
	g++ -g Source2LinuxN.cpp -o Source2LinuxN
	./Source2LinuxN
   Next, enter the parameters as per above.


**4. Running the Analysis**

*   The program will then iterate through each specified chromosome, calculating the various metrics including ABBA and BABA counts and branch lengths in genomic windows of defined size.
*   The output of this will be saved after all chromosomes have been analysed into the specified output `.txt` file.

**5. Setting Run Parameters**

   The program has several modifiable run parameters that are set through user input:
    *   **P1 Population or Region**:  This determines the first group, P1, for comparison.
    *   **P1 Taxon Code**: This is the numerical code of the population or region that you want to compare for P1. Population codes and region codes can be printed to the screen to help users select an appropriate code.
    *   **P2 Population or Region**: This determines the second group for comparison.
    *   **P2 Taxon Code**: This is the numerical code of the population or region that you want to compare for P2. Population codes and region codes can be printed to the screen to help users select an appropriate code.
    *   **Start Chromosome**: This is the chromosome number where the analysis will start.
    *   **End Chromosome**: This is the chromosome number where the analysis will end.
    *   **Window Size**: This is the size of the genomic windows to be analysed, default = 50kb. Radically larger or smaller sizes may have negative impacts on the outputs.
    *   **Spiking Percent**: The proportion of introgressed material used to simulate introgression (entered as a fraction, not a percent!) 

   To modify parameters:
    * The program will print the current settings for P1 and P2, and start and end chromosomes and prompt you to enter a code to modify a parameter.
    *   Enter **`a`** to toggle P1 between population and region.
    *   Enter **`b`** to enter a new taxon code for P1.
    *   Enter **`c`** to toggle P2 between population and region.
    *   Enter **`d`** to enter a new taxon code for P2.
    *   Enter **`e`** to change the starting chromosome.
    *   Enter **`f`** to change the ending chromosome.
    *   Enter **`g`** to change the window size.
    *   Enter **`h`** to change the spiking proportion.
    *   Enter **`Y`** to finish entering parameters and begin the analysis.
    *   Enter **`Z`** to terminate the program.
    * The program will prompt for a new input value for the appropriate parameter.
    * If P1 or P2 are set to region, and the entered code is greater than 6, then it will prompt the user to re-enter the code with a number between 0 and 6.
    *   **Important:** Only consecutive groups of chromosomes can be analyzed.
* The default settings analyze chromosome 22 for the Yoruba and French populations.

**6. Output**

The output is written to a text file, which includes:
    *   D values based on specific subsets of sites (mutating triplet, transitions and transversions, whether the three archaics carry the same of different alleles).
    *   How metrics such as the Neanderthal branch length asymmetry, D and heterozygosity difference vary with the branch length asymmetry to the chimpanzee.
    *   How D varies with the slope of the relationship between heterozygosity and distance from Africa, both overall and based on subsets of sites where either the chimpanzee A or Neanderthal B is the major human allele.
    *   How D varies with African - non-African heterozygosity difference, again, both overall and based on subsets of sites where either the chimpanzee A or Neanderthal B is the major human allele.
    *   The meaning of the figures and table information is included in the comments of the `Generate_output` function.
    *   This code produces: Figures 4, 5, 6, 11 and S3

**7. Key Concepts in the Code**
* The program generates values calculated both in non-overlapping genomic windows and in contiguous blocks of 100,000 informative sites over all chromosomes selected.
* The program calculates the triplet type (e.g., ""ATA_TV"", ""CGC_TS"") at each site, which are used in the output, this uses the humans reference sequence contained in chr[chromosome].fa.
* The program calculates transition and transversion rates.
* The program refers to archaic hominin populations using specific codes (e.g., Vindija, Altai, Denisova), which are used to calculate ABBA and BABA statistics.   
