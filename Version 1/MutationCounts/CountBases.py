#CountBases
# -*- coding: utf-8 -*-
"""
CountBases
------------
This Python script processes genomic data stored in a tab-delimited file, counting the occurrences of alternative alleles (ALT) in each 50,000-base window for each chromosome.

The ALT field represents the alternative nucleotide (A, C, G, T) in a genomic variant. 
The code identifies whether the ALT contains one, two, or three alternative nucleotides (or alleles) separated by commas, 
counts them, and outputs the total number of single, double, and triple alleles for each window in the genome. 
The goal is to analyze the distribution of variant types across different genomic regions.

The sciprt is applicable to both HGDP and dbSNP, see hte lines at the bottom to run each set of files

Name: Eran Elhaik
Date: 5/9/2024
Version: 1.50

Ver 1.10: I changed the window size to 50kb and added comments
Ver 1.50: For the HGDP, i am loading data from a file and use it to filter SNPs
"""

import csv  # Module for reading and writing CSV files
#from collections import defaultdict  # defaultdict is imported but not used

#non overlapping window sizes for the SNP counts
window_size = 50000

# Function to count letters in ALT column
def count_alt_letters(alt):
    """
    Splits the 'alt' string based on commas and counts the number of entries.
    Depending on the number of entries, increments the single, double, or triple count.
    
    Parameters:
    alt (str): The ALT field from the input file containing comma-separated letters.

    Returns:
    tuple: Counts of single, double, and triple occurrences (in that order).
    """
    single_count = 0
    double_count = 0
    triple_count = 0

    # Split the ALT field by commas to get individual variants
    letters = alt.split(',')

    # Count occurrences based on the number of variants
    if len(letters) == 1:
        single_count += 1  # Count of single letters
    elif len(letters) == 2:
        double_count += 1  # Count of two letters
    elif len(letters) == 3:
        triple_count += 1  # Count of three letters
    
    # Return the counts for single, double, and triple occurrences
    return single_count, double_count, triple_count

# Function to read unique chromosome IDs from Chr_Autosome_IDs.txt
def read_chromosome_ids(chromosome_file):
    """
    Reads unique chromosomes from the provided file and returns them as a set.
    
    Parameters:
    chromosome_file (str): Path to the file containing chromosome IDs.

    Returns:
    set: A set containing unique chromosome IDs.
    """
    chromosome_ids = set()
    with open(chromosome_file, 'r') as file:
        for line in file:
            chromosome_ids.add(line.strip())  # Strip any whitespace/newlines and add to set
    return chromosome_ids


# Main function
def main(input_file, chromosome_file, output_file):
    """
    Reads the input tab-delimited file, processes the data in windows of 10,000 bases
    per chromosome, and counts single, double, and triple ALT occurrences per window.
    
    Only processes chromosomes present in the chromosome_file.

    Writes the results to an output file.
    
    Parameters:
    input_file (str): The path to the tab-delimited file to be processed.
    chromosome_file (str): The path to the file containing unique chromosome IDs.
    output_file (str): The path to the output file where the results will be written.
    """
    
    # Read the unique chromosome IDs
    chromosome_ids = read_chromosome_ids(chromosome_file)    
    
    # Initialize a dictionary to track counts for each window
    window_counts = {'single': 0, 'double': 0, 'triple': 0}
    current_window_index = None
    current_chrom = None

    print("Starting to process the input file:", input_file)

    # Open the output file for writing
    with open(output_file, 'w') as out_file:
        out_file.write("Chromosome\tWindow Start\tWindow End\tSingle Count\tDouble Count\tTriple Count\n")  # Header

        # Open the input file and read its contents
        with open(input_file, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            
            # next line
            #next(reader)
    
            # Iterate through each row in the file
            for row in reader:
                # Ensure the row has at least 5 columns (chrom, pos, alt, etc.)
                if len(row) < 5:
                    continue  # Ensure there are enough columns
                
                # Extract the relevant fields from the row
                chrom, pos, _, _, alt = row
                pos = int(pos)
                
                # Only process rows where the chromosome is in our unique chromosome list
                if chrom not in chromosome_ids:
                    print(chrom)
                    print("found a strange chromosome, skipping...")
                    continue
                
                # Calculate the window index: position // 10,000 determines the window (0-based index)
                window_index = (pos - 1) // window_size
    
                # Check if we're processing a new window or chromosome
                if current_window_index is not None and (window_index != current_window_index or chrom != current_chrom):
                    # Write the counts for the previous window (before switching to the new window/chromosome)
                    out_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        current_chrom, 
                        current_window_index * window_size + 1, 
                        (current_window_index + 1) * window_size, 
                        window_counts['single'], 
                        window_counts['double'], 
                        window_counts['triple']
                    ))
    
                    # Reset the counts for the new window
                    window_counts = {'single': 0, 'double': 0, 'triple': 0}
    
                   # If it's a new chromosome, reset the window counts entirely
                if chrom != current_chrom:
                    print(f"Processing new chromosome: {chrom}")
                    current_window_index = None  # Reset the window index
    
                # Update the current window and chromosome
                current_window_index = window_index
                current_chrom = chrom
    
                # Count the letters in ALT
                single_count, double_count, triple_count = count_alt_letters(alt)
    
                # Accumulate counts in the corresponding window
                window_counts['single'] += single_count
                window_counts['double'] += double_count
                window_counts['triple'] += triple_count
                
                # Optional: Uncomment the following lines to print progress for each line processed
                # print("Processed: CHROM={}, POS={}, ALT={}, Window={} - {}, Single Count={}, Double Count={}, Triple Count={}".format(
                #      chrom, pos, alt,
                #      window_index * window_size + 1, (window_index + 1) * window_size,
                #      single_count, double_count, triple_count))
    
                
                # Print progress for each line
    #            print("Processed: CHROM={}, POS={}, ALT={}, "
    #                  "Window={} - {}, Single Count={}, Double Count={}, Triple Count={}".format(
    #                      chrom, pos, alt,
    #                      window_index * window_size + 1, (window_index + 1) * window_size,
    #                      single_count, double_count, triple_count))
    
        # After processing all rows, output the counts for the last window
        if current_window_index is not None:
            out_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                current_chrom, 
                current_window_index * window_size + 1, 
                (current_window_index + 1) * window_size, 
                window_counts['single'], 
                window_counts['double'], 
                window_counts['triple']
            ))
            

# If this script is run as the main program, execute the main function
if __name__ == "__main__":
    
    #The script was applied to dbSNP and HGDP data (unmask each part to run it)
    
    #dbSNP scripts 	
    #input_file = "output_no_indels_common.txt"  # Input your file here
    #chromosome_file = "Chr_Autosome_IDs.txt"  # File containing unique chromosomes
    #output_file = "processed_output_dbSNP.txt"  # Output file to store the results

    #HGDP scripts 	
    input_file = "/lunarc/nobackup/projects/snic2019-34-3/shared_elhaik_lab1/Projects/Amos/HGDP/Allchr_allSNPs.txt"
    chromosome_file = "Chr_Autosome_IDs_HGDP.txt"  # File containing unique chromosomes
    output_file = "processed_output_HGDP.txt"  # Output file to store the results

    main(input_file, chromosome_file, output_file)
