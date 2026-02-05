# Amphipathic_helices

Complete Python implementation of HeliQuest for analyzing amphipathic helices in protein sequences.

## Files

- `heliquest_final.py` - Core analysis module
- `analyze_fasta.py` - Batch processor for FASTA files
- `example_usage.py` - Examples of how to use the module
- `test_sequences.fasta` - Sample FASTA file for testing

## Quick Start

### 1. Analyze a Single Sequence

```python
from heliquest_final import analyze_amphipathic_helices

sequence = "SPLTDSGGYYYPKGELDSLLEDIDSPFSIISSQKGKKS..."
best_sequence, max_moment, net_charge = analyze_amphipathic_helices(sequence)

print(f"Best amphipathic helix: {best_sequence}")
print(f"Hydrophobic moment: {max_moment:.3f}")
print(f"Net charge: {net_charge}")
```

### 2. Analyze a FASTA File

```bash
python analyze_fasta.py your_proteins.fasta
```

This will create a timestamped output directory with:
- Individual analysis for each protein (CSV, plot, summary)
- Cumulative summary CSV
- Statistics report
- HTML report

### 3. Advanced Options

```bash
# Custom window size
python analyze_fasta.py proteins.fasta --window 11

# Custom output directory
python analyze_fasta.py proteins.fasta --output /path/to/results

# View help
python analyze_fasta.py --help
```

## Output Structure

```
fasta_name_TIMESTAMP/
├── cumulative_summary.csv          # All results in one table
├── cumulative_statistics.txt       # Statistical summary
├── report.html                     # Interactive HTML report
└── individual_proteins/
    ├── 001_protein_name/
    │   ├── heliquest_results.csv   # Full window analysis
    │   ├── heliquest_wheel.png     # Helical wheel plot
    │   └── summary.txt             # Text summary
    ├── 002_protein_name/
    │   └── ...
    └── ...
```

## Function Reference

### analyze_amphipathic_helices()

```python
analyze_amphipathic_helices(
    seq,                                    # Protein sequence (str)
    window=18,                              # Window size (int)
    output_dir="/mnt/user-data/outputs",    # Output directory (str)
    save_files=True,                        # Save files? (bool)
    show_plot=False                         # Display plot? (bool)
)
```

**Returns:** `(best_sequence, max_moment, net_charge)`

## Output Files Explained

### cumulative_summary.csv
Contains one row per protein with:
- Protein ID and name
- Sequence length
- Best helix position and sequence
- Hydrophobic moment (μH)
- Hydrophobicity (H)
- Net charge (z)
- Polar/nonpolar frequencies
- Path to individual results

### Individual Protein Folders
Each protein gets its own folder containing:

1. **heliquest_results.csv** - Complete sliding window analysis
   - All possible 18-residue windows
   - Calculated properties for each window

2. **heliquest_wheel.png** - Helical wheel projection
   - Visual representation of the best helix
   - Color-coded residues (blue=K/R, red=D/E, yellow=hydrophobic, etc.)
   - Orange arrow shows hydrophobic moment direction
   - N and C terminus labels

3. **summary.txt** - Text summary of best helix

## Residue Classifications

**Polar (STNQYHDEKR):**
- Polar uncharged: S, T, N, Q
- Aromatic: Y, H
- Charged: D, E, K, R

**Nonpolar (AVLIMFWPCG):**
- Hydrophobic: A, V, L, I, M, F, W
- Special: P, C, G

## Examples

See `example_usage.py` for:
1. Basic single sequence analysis
2. Analysis without saving files
3. Using different window sizes
4. Batch processing multiple sequences

## Algorithm Details

Uses the Fauchère-Pliska hydrophobicity scale and calculates:

- **Mean Hydrophobicity (H)**: Average hydrophobicity of window
- **Hydrophobic Moment (μH)**: Vector sum of hydrophobicities at 100° intervals
- **Net Charge (z)**: (K+R) - (D+E)
- **Residue Frequencies**: Proportion of polar/nonpolar residues

The helical wheel uses α-helix geometry (100° rotation per residue).

## Requirements

```python
numpy
pandas
matplotlib
```

## Validation

The implementation has been validated against HeliQuest output and matches perfectly on all calculated values (hydrophobicity, hydrophobic moment, net charge, and residue frequencies).

## Citation

If you use this tool, please cite the original HeliQuest:
Gautier R, et al. (2008) HELIQUEST: a web server to screen sequences with specific alpha-helical properties. Bioinformatics 24(18):2101-2102.

## Tips

- Use window size 18 for typical amphipathic helices
- Use window size 11 for shorter helical segments
- Check the HTML report for an interactive overview
- The hydrophobic moment arrow in the wheel plot shows the amphipathic "direction"
- Longer arrows = stronger amphipathicity
- Sequences shorter than the window size are automatically skipped

## Troubleshooting

**"Sequence too short" errors**: Use a smaller window size or filter sequences in your FASTA file

**Import errors**: Make sure `heliquest_final.py` is in the same directory or in your Python path

**No output files**: Check that `save_files=True` and the output directory is writable