# PlotWindowCounts
# -*- coding: utf-8 -*-
"""
PlotWindowCounts Script
=======================

This script processes and visualizes genomic mutation data, focusing on the percentage distribution of mutation types (Single, Double, and Triple Counts) across genomic windows for multiple chromosomes. Key functionalities include:

1. **Data Input**: Reads a tab-delimited file containing genomic data into a pandas DataFrame.
2. **Window Aggregation**: Groups smaller genomic windows into larger aggregated windows, summing mutation counts for Single, Double, and Triple mutations.
3. **Multi-Chromosome Visualization**: Generates a grid of stacked bar charts, each representing the mutation distribution for a specific chromosome, with the data dynamically aggregated based on window size.
4. **Output**: Saves a consolidated plot as a high-resolution PNG file.

Notes:
- Processes and visualizes all chromosomes in a single execution.
- Automatically adjusts the number of rows and columns for the plot grid.
- Removes unused subplots for cleaner visualization.

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

# Constants
windows_in_screen = 200

# Read the file into a pandas dataframe
file_path = "processed_output_dbSNP_50k.txt"  # Assuming this file exists
file_path = "processed_output_all_chrs_HGDP_50k.txt"  # Assuming this file exists
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
        # Select the current group of windows
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

# Function to plot percentage of mutations per window for all chromosomes
def plot_all_chromosomes(data):
    # Create a figure and axis with 2 columns
    num_chromosomes = data['Chromosome'].nunique()
#    num_rows = (num_chromosomes + 1) // 2  # Calculate rows needed for 2 plots per row
#    fig, axes = plt.subplots(num_rows, 2, figsize=(15, 5 * num_rows))

    num_cols = 4  # Set 4 columns
    num_rows = (num_chromosomes + num_cols - 1) // num_cols  # Calculate rows needed for 4 plots per row
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 5 * num_rows))  # Adjust figure size accordingly

    # Flatten the axes array for easier indexing
    axes = axes.flatten()

    # Iterate through each chromosome and plot
    for idx, (chromosome, group) in enumerate(data.groupby('Chromosome')):
        print(f"Now analyzing chromosome {chromosome} with {len(group)} rows")
        
        # Aggregate the data every X windows
        window_size = round(len(group) / windows_in_screen)
        aggregated_data = aggregate_windows(group, window_size)

        # Create window labels for the X-axis
        window_labels = aggregated_data.apply(lambda row: f"{row['Window Start']}-{row['Window End']}", axis=1)
        num_windows = len(window_labels)
        x = np.arange(num_windows)

        # Calculate total counts and percentages
        aggregated_data['Total Count'] = aggregated_data[['Single Count', 'Double Count', 'Triple Count']].sum(axis=1)
        aggregated_data['Single %'] = aggregated_data['Single Count'] / aggregated_data['Total Count'] * 100
        aggregated_data['Double %'] = aggregated_data['Double Count'] / aggregated_data['Total Count'] * 100
        aggregated_data['Triple %'] = aggregated_data['Triple Count'] / aggregated_data['Total Count'] * 100

        # Plot stacked bars
        axes[idx].bar(x, aggregated_data['Single %'], label='Single Count', color='blue', alpha=1)
        axes[idx].bar(x, aggregated_data['Double %'], label='Double Count', color='orange', alpha=1, bottom=aggregated_data['Single %'])
        axes[idx].bar(x, aggregated_data['Triple %'], label='Triple Count', color='black', alpha=1, bottom=aggregated_data['Single %'] + aggregated_data['Double %'])

        # Set labels and title
        #axes[idx].set_ylabel('Percentage of Mutations (%)', fontsize=12)
        #axes[idx].set_title(f'Chromosome {chromosome}', fontsize=14)

        # Set x-ticks to show every 5th label
        '''
        if num_windows > 10:
            axes[idx].set_xticks(np.arange(0, num_windows, 10))  # Set ticks at every 5th index
            axes[idx].set_xticklabels(window_labels[::10], rotation=45, ha="right")  # Set labels for every 5th tick
'''
        # Inside the plot_all_chromosomes function, after plotting the bars
        axes[idx].set_xlim(0, num_windows)  # Set x-limits for the current subplot

        # Set y-axis limits to 0-100%
        axes[idx].set_ylim(0, 100)

        #Hide x axis numbers
        axes[idx].set_xticks([])  # Removes the ticks and labels on the x-axis

        # Add legend
        #axes[idx].legend()

    # Remove any empty subplots
    for i in range(idx + 1, len(axes)):
        fig.delaxes(axes[i])

    # Adjust layout to fit plots closely
    plt.tight_layout()
    plt.savefig('all_chromosomes_plot.png', dpi=600)
    plt.show()

# Call the function to plot all chromosomes
plot_all_chromosomes(data)

print("Program PlotWindowCounts completed.")
