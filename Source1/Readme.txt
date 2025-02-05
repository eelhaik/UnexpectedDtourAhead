**Instructions for Using Source.cpp**

**1. Introduction**

This program calculates aspects of the introgression statistic D, broken down into transitions and transversions (TS/TV), which 3-base motif (triplet) has mutated and the frequencies of the alleles in different populations. 
The output generated is a tab delimited file that is best viewed by opening in Excel or equivalent and start with two blocks of overall state counts.
The first block contains counts of different base conformations partitioned by TS/TV and triplet; the second contains partitioned ABBA and BABA counts for D(P1,P2,archaic,chimpanzee). 
In the second block, ABBA and BABA counts are partitioned by TS/TV and triplet (rows) and by which allele is the major allele in humans and which alleles fixed / polymorphic in each human population P1 and P2.
Thus, opening the output file in Excel, column D contains BABA counts for D(P1,P2,Vindija,chimpanzee) where the human major allele is A and P/A indicates polymorphism in P1 and P2 fixed for A.
Several columns contain either all or mainly zeros reflecting the fact that some scenarios make either ABBAs or BABAs impossible or extremely unlikely (when P2 is fixed for A, ABBA cannot occur).
Below these full outputs are a series of Tables corresponding to specific figures in the paper. Explanations for each are in the text of the paper.
Note, analysing just the smallest chromosome 22 results in zeros in some cells for the rarestbase combinations.

**2. Prerequisites**

*   The program is designed to be compiled and run in a C++ environment.
*   The following input files should be located in the same directory:
    *   HGDP VCF files (e.g., `hgdp.v0.5.archaics.chr22.vcf`)
    *   Reference FASTA files (e.g., `chr22.fa`)
    *   Population information file (`inpopsHGDP.txt`)
*   Filenames must  be exact and are case-sensitive.

**3. Initial Setup**

*   The program uses a **graphical interface** to select the directory with the input files. The program will prompt you to select a folder using a dialog box.
*   After selecting the folder, the program will load the population data from `inpopsHGDP.txt` and the reference sequence data from the `chr[chromosome number].fa` files.
*   The program will then prompt you for an **output file name**. The output file will be a `.txt` file and is saved into the selected directory.

* In Linux, there is no graphic interface. To compile the code do:
	g++ -g SourceLinux.cpp -o SourceLinux
	./SourceLinux
   Next, enter the parameters as per above.

**4. Running the Analysis**

*   The program will then iterate through each specified chromosome, calculates the various metrics including ABBA and BABA counts.
*   The output of this will be saved after all chromosomes have been anlysed into the specified output `.txt` file.

**5. Setting Run Parameters**

   The program has several modifiable run parameters that are set through user input:
    *   **P1 Population or Region**:  This determines the first group, P1, for comparison.
    *   **P1 Taxon Code**: This is the numerical code of the population or region that you want to compare for P1. Population codes and region codes can be printed to the screen to help users select an appropriate code.
    *   **P2 Population or Region**: This determines the second group for comparison.
    *   **P2 Taxon Code**: This is the numerical code of the population or region that you want to compare for P2. Population codes and region codes can be printed to the screen to help users select an appropriate code.
    *   **Start Chromosome**: This is the chromosome number where the analysis will start.
    *   **End Chromosome**: This is the chromosome number where the analysis will end.

	Parameter values for each run are given at the top of the output file.

   To modify parameters:
    * The program will print the current settings for P1 and P2, and start and end chromosomes and prompt you to enter a code to modify a parameter.
    *   Enter **`a`** to toggle P1 between population and region.
    *   Enter **`b`** to enter a new taxon code for P1.
    *   Enter **`c`** to toggle P2 between population and region.
    *   Enter **`d`** to enter a new taxon code for P2.
    *   Enter **`e`** to change the starting chromosome.
    *   Enter **`f`** to change the ending chromosome.
    *   Enter **`Y`** to finish entering parameters and begin the analysis.
    *   Enter **`Z`** to terminate the program.
    * The program will prompt for a new input value for the appropriate parameter.
    * If P1 or P2 are set to region, and the entered code is greater than 6, then it will prompt the user to re-enter the code with a number between 0 and 6.
    *   **Important:** Only consecutive groups of chromosomes can be analyzed.
* The default settings analyze chromosome 22 for the Yoruba and French populations.

**6. Output**

The output is written to a text file, which includes:
    *   Counts of sites by class, including triplet types and counts for various categories of differences (Vindija, Altai, Denisovan, Neanderthals, and humans).
    *   ABBA and BABA counts.
    *   Tables summarizing counts and ratios of different types of sites.
    *   Figures (as text) displaying proportions and differences in ABBA/BABA counts and other metrics.
    *   Tables with D* statistics, comparing different archaic groups.
    *   The meaning of the figures and table information is included in the comments of the `Generate_output` function.
    *   This code produces: Tables 2 and 3 and Figures 7, 8 and 9, 10 and S2

**7. Key Concepts in the Code**
* The program generates pooled values calculated over all chromosomes selected.
* The program calculates the triplet type (e.g., "ATA_TV", "CGC_TS") at each site, which are used in the output, this uses the humans reference sequence contained in chr[chromosome].fa.
* The program calculates transition and transversion rates.
* The program refers to archaic hominin populations using specific codes (e.g., Vindija, Altai, Denisova), which are used to calculate ABBA and BABA statistics.

