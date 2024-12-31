#PlotWindowCounts
# -*- coding: utf-8 -*-

"""
PlotWindowCounts Script
=======================

This script analyzes genomic data from a tab-delimited file, aggregating mutation counts across genomic windows and visualizing the percentage distribution of different mutation types (Single, Double, and Triple Counts) for each chromosome. It performs the following tasks:

1. **Data Input**: Reads genomic data from a specified input file into a pandas DataFrame.
2. **Aggregation**: Groups genomic windows into larger aggregated windows, summing mutation counts for each category.
3. **Visualization**: Generates stacked bar charts displaying the percentage of Single, Double, and Triple Counts for each aggregated window.
4. **Chromosome-wise Analysis**: Processes each chromosome separately, dynamically calculating window sizes for aggregation.

Output includes high-resolution plots for each chromosome, saved as PNG files.

Dependencies:
- `matplotlib`
- `pandas`
- `numpy`

Run it: Change the file_path to get the results for dbSNP and HGDP datasets

Author: Eran Elhaik
Created on: Thu Sep 5, 2024
Ver: 1.00
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys

#constants
windows_in_screen = 200

# Read the file into a pandas dataframe
file_path = "processed_output_dbSNP_50k.txt"  # Assuming this file exists
#file_path = "processed_output_all_chrs_HGDP_50k.txt"  # Assuming this file exists
data = pd.read_csv(file_path, sep="\t")

print("Start program PlotWindowCounts.")

# Function to aggregate every X windows and calculate their sum
def aggregate_windows(data, window_size):
    # Number of windows
    num_windows = len(data)

    # Calculate how many aggregated windows we will have
    num_aggregates = (num_windows + window_size - 1) // window_size  # Round up

    # List to hold the aggregated data
    aggregated_data = []

    for i in range(num_aggregates):
        # Select the current group of 20 windows
        start_idx = i * window_size
        end_idx = min((i + 1) * window_size, num_windows)

        # Get the subset of data
        window_subset = data.iloc[start_idx:end_idx]

        # Calculate the sum of the counts for this subset
        sum_single_count = window_subset['Single Count'].sum()
        sum_double_count = window_subset['Double Count'].sum()
        sum_triple_count = window_subset['Triple Count'].sum()

        # Define start and end window labels
        start_window = window_subset['Window Start'].min()
        end_window = window_subset['Window End'].max()

        # Append the aggregated result
        aggregated_data.append({
            'Window Start': start_window,
            'Window End': end_window,
            'Single Count': sum_single_count,
            'Double Count': sum_double_count,
            'Triple Count': sum_triple_count
        })

    # Convert to a pandas DataFrame
    return pd.DataFrame(aggregated_data)


# Function to plot grouped bar graphs for each window
def plot_counts(data, chromosome):
    
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(15, 6))

    # Create window labels for the X-axis
    window_labels = data.apply(lambda row: f"{row['Window Start']}-{row['Window End']}", axis=1)

    # Calculate the number of windows
    num_windows = len(window_labels)

    # Generate positions for the bars
    x = np.arange(num_windows)

    # Calculate total counts and percentages
    data['Total Count'] = data[['Single Count', 'Double Count', 'Triple Count']].sum(axis=1)
    data['Single %'] = data['Single Count'] / data['Total Count'] * 100
    data['Double %'] = data['Double Count'] / data['Total Count'] * 100
    data['Triple %'] = data['Triple Count'] / data['Total Count'] * 100

    # Plot stacked bars
    ax.bar(x, data['Single %'], label='Single Count', color='blue', alpha=1)
    ax.bar(x, data['Double %'], label='Double Count', color='orange', alpha=1, bottom=data['Single %'])
    ax.bar(x, data['Triple %'], label='Triple Count', color='black', alpha=1, bottom=data['Single %'] + data['Double %'])

    # Set labels and title
    #ax.set_xlabel('Windows (Start-End)', fontsize=12)
    #ax.set_ylabel('Percentage of Mutations (%)', fontsize=12)

    # Set y-axis limits to 0-100%
    ax.set_ylim(0, 100)
    
    # Hide the x-axis
    #ax.xaxis.set_visible(False)
    
    # Set X-axis ticks and labels
    #ax.set_xticks(x)
    #ax.set_xticklabels(window_labels, rotation=45, ha="right")
    
    '''
    # Set x-ticks to show every 10th label
    if num_windows > 5:
        ax.set_xticks(np.arange(0, num_windows, 5))  # Set ticks at every 10th index
        ax.set_xticklabels(window_labels[::5], rotation=45, ha="right")  # Set labels for every 10th tick
'''
    
    # Set limits to remove gaps between bars and the y-axis
    ax.set_xlim(0, num_windows)

    #Hide x axis numbers
    ax.get_xaxis().set_visible(False)
    
    # Add legend
    #ax.legend()

    # Save the figure at 600 dpi before showing
    print("Saving figure for chromosome:", chromosome, "...")
    plt.tight_layout()
    plt.savefig(f'{chromosome}_figure.png', dpi=600)
    plt.show()

# Group the data by Chromosome and analyze each separately
for chromosome, group in data.groupby('Chromosome'):
    print(f"Now analyzing chromosome {chromosome} with {len(group)} rows")
    
    # Aggregate the data every X windows
    window_size = round(len(group)/windows_in_screen)
    aggregated_data = aggregate_windows(group, window_size)

    # Save the aggregated data to a file named after the chromosome
    #aggregated_data.to_csv(f'{chromosome}_aggregated.csv', index=False, sep="\t")

    # Call the function to plot the data
    plot_counts(aggregated_data, chromosome)
    
    sys.exit()


print("Program PlotWindowCounts completed.")


