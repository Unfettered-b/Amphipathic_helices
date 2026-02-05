#!/usr/bin/env python3
"""
FASTA Amphipathic Helix Analyzer

Process a FASTA file and analyze all sequences for amphipathic helices.
Generates individual reports for each protein and a cumulative summary.

Usage:
    python analyze_fasta.py input.fasta
    python analyze_fasta.py input.fasta --window 18 --output /custom/path
"""

import sys
import os
import argparse
from datetime import datetime
from pathlib import Path
import pandas as pd

# Import the heliquest analyzer
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from heliquest_final import analyze_amphipathic_helices, analyze_sequence


def parse_fasta(fasta_file):
    """
    Parse a FASTA file and return a list of (header, sequence) tuples.
    
    Parameters:
    -----------
    fasta_file : str
        Path to FASTA file
        
    Returns:
    --------
    list : [(header, sequence), ...]
    """
    sequences = []
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                # Add to current sequence (remove any whitespace)
                current_seq.append(line.replace(' ', '').replace('\t', ''))
        
        # Don't forget the last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))
    
    return sequences


def clean_filename(name):
    """Convert protein name to valid filename"""
    # Remove or replace problematic characters
    name = name.split()[0]  # Take first word
    name = name.replace('|', '_').replace('/', '_').replace('\\', '_')
    name = name.replace('>', '').replace('<', '').replace(':', '_')
    name = name.replace('*', '').replace('?', '').replace('"', '')
    return name[:50]  # Limit length


def analyze_fasta_file(fasta_file, window=18, output_base_dir=None):
    """
    Analyze all sequences in a FASTA file for amphipathic helices.
    
    Parameters:
    -----------
    fasta_file : str
        Path to FASTA file
    window : int
        Window size for helix analysis (default: 18)
    output_base_dir : str or None
        Base directory for outputs. If None, creates in current directory.
        
    Returns:
    --------
    pd.DataFrame : Summary of all proteins analyzed
    """
    
    # Create output directory with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fasta_basename = Path(fasta_file).stem
    
    if output_base_dir is None:
        output_base_dir = "."
    
    output_dir = os.path.join(output_base_dir, f"{fasta_basename}_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 80)
    print(f"FASTA AMPHIPATHIC HELIX ANALYZER")
    print("=" * 80)
    print(f"Input file:     {fasta_file}")
    print(f"Output dir:     {output_dir}")
    print(f"Window size:    {window}")
    print(f"Timestamp:      {timestamp}")
    print("=" * 80)
    
    # Parse FASTA file
    sequences = parse_fasta(fasta_file)
    print(f"\nFound {len(sequences)} sequences in FASTA file")
    
    # Filter sequences >= window size
    valid_sequences = [(h, s) for h, s in sequences if len(s) >= window]
    skipped = len(sequences) - len(valid_sequences)
    
    if skipped > 0:
        print(f"Skipped {skipped} sequences shorter than {window} residues")
    print(f"Analyzing {len(valid_sequences)} sequences...\n")
    
    # Create subdirectories
    individual_dir = os.path.join(output_dir, "individual_proteins")
    os.makedirs(individual_dir, exist_ok=True)
    
    # Store results for cumulative summary
    cumulative_results = []
    
    # Process each sequence
    for idx, (header, sequence) in enumerate(valid_sequences, 1):
        protein_name = clean_filename(header)
        print(f"[{idx}/{len(valid_sequences)}] Processing: {header[:60]}...")
        
        try:
            # Create protein-specific subdirectory
            protein_dir = os.path.join(individual_dir, f"{idx:03d}_{protein_name}")
            os.makedirs(protein_dir, exist_ok=True)
            
            # Analyze the sequence
            best_seq, max_moment, net_charge = analyze_amphipathic_helices(
                sequence,
                window=window,
                output_dir=protein_dir,
                save_files=True,
                show_plot=False
            )
            
            # Get full analysis for this protein
            df = analyze_sequence(sequence)
            best_idx = df["HydrophobicMoment"].idxmax()
            best_row = df.loc[best_idx]
            
            # Save detailed results
            result_info = {
                'Protein_ID': idx,
                'Protein_Name': header,
                'Sequence_Length': len(sequence),
                'Best_Helix_Position': f"{best_row['Start']}-{best_row['End']}",
                'Best_Helix_Sequence': best_seq,
                'Hydrophobic_Moment': max_moment,
                'Hydrophobicity': float(best_row['Hydrophobicity']),
                'Net_Charge': net_charge,
                'Freq_Polar': float(best_row['FreqPolar']),
                'Freq_NonPolar': float(best_row['FreqNonPolar']),
                'Output_Dir': protein_dir
            }
            cumulative_results.append(result_info)
            
            # Save protein-specific summary
            with open(os.path.join(protein_dir, "summary.txt"), 'w') as f:
                f.write(f"Protein: {header}\n")
                f.write(f"Sequence Length: {len(sequence)} residues\n")
                f.write(f"Window Size: {window}\n\n")
                f.write(f"MOST AMPHIPATHIC HELIX:\n")
                f.write(f"  Position: {best_row['Start']}-{best_row['End']}\n")
                f.write(f"  Sequence: {best_seq}\n")
                f.write(f"  Hydrophobic Moment: {max_moment:.3f}\n")
                f.write(f"  Hydrophobicity: {best_row['Hydrophobicity']:.3f}\n")
                f.write(f"  Net Charge: {net_charge}\n")
                f.write(f"  Freq Polar: {best_row['FreqPolar']:.3f}\n")
                f.write(f"  Freq NonPolar: {best_row['FreqNonPolar']:.3f}\n")
            
            print(f"    → Best helix: {best_seq} (μH={max_moment:.3f}, z={net_charge})")
            
        except Exception as e:
            print(f"    ✗ ERROR: {str(e)}")
            cumulative_results.append({
                'Protein_ID': idx,
                'Protein_Name': header,
                'Sequence_Length': len(sequence),
                'Best_Helix_Position': 'ERROR',
                'Best_Helix_Sequence': 'ERROR',
                'Hydrophobic_Moment': None,
                'Hydrophobicity': None,
                'Net_Charge': None,
                'Freq_Polar': None,
                'Freq_NonPolar': None,
                'Output_Dir': 'ERROR: ' + str(e)
            })
    
    # Create cumulative summary
    print("\n" + "=" * 80)
    print("Creating cumulative summary...")
    
    summary_df = pd.DataFrame(cumulative_results)
    
    # Save cumulative CSV
    cumulative_csv = os.path.join(output_dir, "cumulative_summary.csv")
    summary_df.to_csv(cumulative_csv, index=False)
    print(f"✓ Cumulative CSV saved: {cumulative_csv}")
    
    # Save cumulative statistics
    stats_file = os.path.join(output_dir, "cumulative_statistics.txt")
    with open(stats_file, 'w') as f:
        f.write("CUMULATIVE STATISTICS\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Total proteins analyzed: {len(valid_sequences)}\n")
        f.write(f"Window size: {window}\n")
        f.write(f"Timestamp: {timestamp}\n\n")
        
        # Calculate statistics (excluding errors)
        valid_results = summary_df[summary_df['Hydrophobic_Moment'].notna()]
        
        if len(valid_results) > 0:
            f.write(f"HYDROPHOBIC MOMENT STATISTICS:\n")
            f.write(f"  Mean:   {valid_results['Hydrophobic_Moment'].mean():.3f}\n")
            f.write(f"  Median: {valid_results['Hydrophobic_Moment'].median():.3f}\n")
            f.write(f"  Std:    {valid_results['Hydrophobic_Moment'].std():.3f}\n")
            f.write(f"  Min:    {valid_results['Hydrophobic_Moment'].min():.3f}\n")
            f.write(f"  Max:    {valid_results['Hydrophobic_Moment'].max():.3f}\n\n")
            
            f.write(f"NET CHARGE STATISTICS:\n")
            f.write(f"  Mean:   {valid_results['Net_Charge'].mean():.2f}\n")
            f.write(f"  Median: {valid_results['Net_Charge'].median():.1f}\n")
            f.write(f"  Min:    {valid_results['Net_Charge'].min()}\n")
            f.write(f"  Max:    {valid_results['Net_Charge'].max()}\n\n")
            
            # Top 10 most amphipathic helices
            f.write(f"TOP 10 MOST AMPHIPATHIC HELICES:\n")
            f.write("-" * 80 + "\n")
            top10 = valid_results.nlargest(10, 'Hydrophobic_Moment')
            for idx, row in top10.iterrows():
                f.write(f"{row['Protein_ID']:3d}. {row['Protein_Name'][:50]}\n")
                f.write(f"     Sequence: {row['Best_Helix_Sequence']}\n")
                f.write(f"     μH={row['Hydrophobic_Moment']:.3f}, z={row['Net_Charge']}\n\n")
    
    print(f"✓ Statistics saved: {stats_file}")
    
    # Create HTML report (optional, nice to have)
    html_file = os.path.join(output_dir, "report.html")
    create_html_report(summary_df, html_file, fasta_file, window, timestamp)
    print(f"✓ HTML report saved: {html_file}")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"Results saved in: {output_dir}")
    print(f"  - Individual protein analyses: {individual_dir}")
    print(f"  - Cumulative summary: {cumulative_csv}")
    print(f"  - Statistics: {stats_file}")
    print(f"  - HTML report: {html_file}")
    print("=" * 80)
    
    return summary_df


def create_html_report(summary_df, html_file, fasta_file, window, timestamp):
    """Create an HTML report of the analysis"""
    
    valid_results = summary_df[summary_df['Hydrophobic_Moment'].notna()]
    
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Amphipathic Helix Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }}
        .header {{ background-color: #2c3e50; color: white; padding: 20px; border-radius: 5px; }}
        .section {{ background-color: white; margin: 20px 0; padding: 20px; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        table {{ border-collapse: collapse; width: 100%; }}
        th {{ background-color: #34495e; color: white; padding: 10px; text-align: left; }}
        td {{ padding: 8px; border-bottom: 1px solid #ddd; }}
        tr:hover {{ background-color: #f5f5f5; }}
        .stat {{ display: inline-block; margin: 10px 20px; }}
        .stat-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
        .stat-label {{ font-size: 12px; color: #7f8c8d; }}
        .sequence {{ font-family: 'Courier New', monospace; background-color: #ecf0f1; padding: 2px 5px; border-radius: 3px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Amphipathic Helix Analysis Report</h1>
        <p><strong>FASTA File:</strong> {fasta_file}</p>
        <p><strong>Analysis Date:</strong> {timestamp}</p>
        <p><strong>Window Size:</strong> {window} residues</p>
    </div>
    
    <div class="section">
        <h2>Summary Statistics</h2>
        <div class="stat">
            <div class="stat-value">{len(summary_df)}</div>
            <div class="stat-label">Proteins Analyzed</div>
        </div>
        <div class="stat">
            <div class="stat-value">{valid_results['Hydrophobic_Moment'].mean():.3f}</div>
            <div class="stat-label">Mean μH</div>
        </div>
        <div class="stat">
            <div class="stat-value">{valid_results['Hydrophobic_Moment'].max():.3f}</div>
            <div class="stat-label">Max μH</div>
        </div>
        <div class="stat">
            <div class="stat-value">{valid_results['Net_Charge'].mean():.1f}</div>
            <div class="stat-label">Mean Charge</div>
        </div>
    </div>
    
    <div class="section">
        <h2>Top 10 Most Amphipathic Helices</h2>
        <table>
            <tr>
                <th>Rank</th>
                <th>Protein</th>
                <th>Sequence</th>
                <th>μH</th>
                <th>Charge</th>
                <th>Position</th>
            </tr>
"""
    
    top10 = valid_results.nlargest(10, 'Hydrophobic_Moment')
    for rank, (idx, row) in enumerate(top10.iterrows(), 1):
        html += f"""
            <tr>
                <td>{rank}</td>
                <td>{row['Protein_Name'][:50]}</td>
                <td class="sequence">{row['Best_Helix_Sequence']}</td>
                <td>{row['Hydrophobic_Moment']:.3f}</td>
                <td>{row['Net_Charge']}</td>
                <td>{row['Best_Helix_Position']}</td>
            </tr>
"""
    
    html += """
        </table>
    </div>
    
    <div class="section">
        <h2>All Proteins</h2>
        <table>
            <tr>
                <th>ID</th>
                <th>Protein</th>
                <th>Length</th>
                <th>Best Sequence</th>
                <th>μH</th>
                <th>Charge</th>
            </tr>
"""
    
    for idx, row in summary_df.iterrows():
        html += f"""
            <tr>
                <td>{row['Protein_ID']}</td>
                <td>{row['Protein_Name'][:50]}</td>
                <td>{row['Sequence_Length']}</td>
                <td class="sequence">{row['Best_Helix_Sequence']}</td>
                <td>{row['Hydrophobic_Moment'] if pd.notna(row['Hydrophobic_Moment']) else 'N/A'}</td>
                <td>{row['Net_Charge'] if pd.notna(row['Net_Charge']) else 'N/A'}</td>
            </tr>
"""
    
    html += """
        </table>
    </div>
</body>
</html>
"""
    
    with open(html_file, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(
        description='Analyze FASTA file for amphipathic helices',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python analyze_fasta.py proteins.fasta
    python analyze_fasta.py proteins.fasta --window 11
    python analyze_fasta.py proteins.fasta --output /path/to/results
        """
    )
    
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--window', '-w', type=int, default=18,
                        help='Window size for helix analysis (default: 18)')
    parser.add_argument('--output', '-o', default='.',
                        help='Base output directory (default: current directory)')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.fasta_file):
        print(f"ERROR: File not found: {args.fasta_file}")
        sys.exit(1)
    
    # Run analysis
    try:
        analyze_fasta_file(args.fasta_file, window=args.window, output_base_dir=args.output)
    except Exception as e:
        print(f"\nERROR: Analysis failed: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()