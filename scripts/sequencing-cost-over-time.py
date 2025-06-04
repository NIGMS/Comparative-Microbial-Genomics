#! /usr/bin/env python3

### 1.2.2 Generate plot showing the cost of sequencing genomes over time

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Load the CSV file
data_table = "data/Sequencing_Cost_Data_Table_May2022.csv"
df = pd.read_csv(data_table)

# Convert the 'Date' column to datetime format using YYYY-MM
df['Date'] = pd.to_datetime(df['Date'], format='%Y-%m')

# Clean up 'Cost per Genome' column by removing $ and commas, and converting to numeric
df['Cost per Genome'] = df['Cost per Genome'].replace({'\$': '', ',': ''}, regex=True).astype(float)

# Plot the data
plt.figure(figsize=(10, 6))

# X-axis: Date, Y-axis: Cost per Genome in log scale
plt.plot(df['Date'], df['Cost per Genome'], marker='o', linestyle='-', color='b')

# Set log scale for y-axis
plt.yscale('log')

# Set labels and title
plt.xlabel('Year', fontsize=14, fontweight='bold')
plt.ylabel('Cost per Genome (Log scale)', fontsize=12, fontweight='bold')
plt.title('Cost per Genome Over Time', fontsize=14, fontweight='bold')

# Format the x-axis for better date presentation
import matplotlib.dates as mdates

# Set major ticks to display years at specific intervals
plt.gca().xaxis.set_major_locator(mdates.YearLocator(3))  # Adjust the number here for the interval
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

# Rotate and format x-ticks
plt.xticks(rotation=45, fontsize=8, fontweight='bold')


# Set the y-axis limits from $100 to $100,000
# plt.ylim(100, 0000)
plt.yticks(fontsize=12, fontweight='bold')

# Display the plot
plt.tight_layout()
plt.savefig("cost_per_genome_plot.png", dpi=300)
plt.show()