import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

WINDOW = 18
ANGLE_RAD = math.radians(100.0)

# Fauchère–Pliska hydrophobicity scale
HYDRO = {
    'A': 0.31,  'R': -1.01, 'N': -0.60, 'D': -0.77,
    'C': 1.54,  'Q': -0.22, 'E': -0.64, 'G': 0.00,
    'H': 0.13,  'I': 1.80,  'L': 1.70,  'K': -0.99,
    'M': 1.23,  'F': 1.79,  'P': 0.72,  'S': -0.04,
    'T': 0.26,  'W': 2.25,  'Y': 0.96,  'V': 1.22
}

# HeliQuest's residue classifications (discovered by reverse engineering)
POLAR = set("STNQYHDEKR")      # Polar uncharged + Y, H + charged residues
NONPOLAR = set("AVLIMFWPCG")   # Hydrophobic + G

def mean_hydrophobicity(seg):
    return sum(HYDRO[a] for a in seg) / len(seg)

def hydrophobic_moment(seg):
    """Calculate hydrophobic moment vector"""
    hx, hy = 0.0, 0.0
    theta = 0.0
    for aa in seg:
        h = HYDRO[aa]
        hx += h * math.cos(theta)
        hy += h * math.sin(theta)
        theta += ANGLE_RAD
    mu = math.sqrt(hx**2 + hy**2) / len(seg)
    return mu, hx / len(seg), hy / len(seg)

def net_charge(seg):
    return sum(1 if a in "KR" else -1 if a in "DE" else 0 for a in seg)

def residue_freq(seg):
    n = len(seg)
    polar_count = sum(1 for a in seg if a in POLAR)
    nonpolar_count = sum(1 for a in seg if a in NONPOLAR)
    return polar_count / n, nonpolar_count / n

def analyze_sequence(seq):
    rows = []
    for i in range(len(seq) - WINDOW + 1):
        seg = seq[i:i+WINDOW]
        H = mean_hydrophobicity(seg)
        muH, vx, vy = hydrophobic_moment(seg)
        z = net_charge(seg)
        fpol, fnon = residue_freq(seg)
        
        rows.append({
            "Helix": i+1,
            "Start": i+1,
            "End": i+WINDOW,
            "Sequence": seg,
            "Hydrophobicity": round(H, 3),
            "HydrophobicMoment": round(muH, 3),
            "NetCharge": z,
            "FreqPolar": round(fpol, 3),
            "FreqNonPolar": round(fnon, 3),
            "vx": vx,
            "vy": vy
        })
    return pd.DataFrame(rows)

def res_color(a):
    """Color scheme matching HeliQuest"""
    if a in "KR": return "blue"
    if a in "DE": return "red"
    if a in "STNQYH": return "purple"
    if a == "G": return "pink"
    if a == "P": return "limegreen"
    return "yellow"

def plot_wheel(seg, vx, vy, title="Helical Wheel"):
    """Plot helical wheel matching HeliQuest layout"""
    n = len(seg)
    radius = 1.0
    
    # Start at top, go clockwise (matching HeliQuest)
    offset = math.pi / 2
    
    angles = []
    x = []
    y = []
    for i in range(n):
        angle = offset - i * ANGLE_RAD  # Subtract for clockwise
        angles.append(angle)
        x.append(radius * math.cos(angle))
        y.append(radius * math.sin(angle))
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Draw circle
    circle = plt.Circle((0, 0), radius, fill=False, color='black', linewidth=2)
    ax.add_patch(circle)
    
    # Draw connections (i to i+3)
    for i in range(n):
        j = (i + 3) % n
        ax.plot([x[i], x[j]], [y[i], y[j]], color='gray', lw=0.5, alpha=0.4, zorder=1)
    
    # Draw residues
    for i, aa in enumerate(seg):
        color = res_color(aa)
        circle = plt.Circle((x[i], y[i]), 0.12, color=color, ec='black', lw=2, zorder=3)
        ax.add_patch(circle)
        
        txt_color = 'white' if color in ['blue', 'red', 'purple'] else 'black'
        ax.text(x[i], y[i], aa, ha='center', va='center', 
                fontsize=12, weight='bold', color=txt_color, zorder=4)
    
    # Draw hydrophobic moment vector
    # The vector is calculated in a coordinate system where residue i is at angle i*100°
    # The wheel displays residue i at angle (offset - i*100°)
    # So we need to apply the transformation: wheel_angle = offset - calc_angle
    # For the vector: if it points at angle α in calc coords, it points at (offset - α) in wheel coords
    # But wait - that's not right for a vector, we need to think about this differently
    
    # Actually, we're rotating the entire coordinate system by (offset) and flipping direction (-)
    # This is equivalent to rotating by (offset) and then reflecting across the y-axis
    # OR just rotating by -offset backwards
    # The correct transformation is: rotate by (90° - 0°) = 90° and flip x
    # Which means: new_angle = offset - old_angle, OR apply rotation matrix
    
    vec_angle = math.atan2(vy, vx)
    vec_magnitude = math.sqrt(vx**2 + vy**2)
    
    # The transformation is: wheel_angle = offset - calc_angle
    # This is equivalent to: reflect across y-axis, then rotate by offset
    # OR: rotate by offset, then reflect across vertical line
    # For a vector: (x,y) -> (-x, y) then rotate by offset
    # Actually simpler: new_angle = offset - old_angle
    rotated_vec_angle = offset - vec_angle
    
    arrow_scale = 1.5
    arrow_x = arrow_scale * vec_magnitude * math.cos(rotated_vec_angle)
    arrow_y = arrow_scale * vec_magnitude * math.sin(rotated_vec_angle)
    
    ax.arrow(0, 0, arrow_x, arrow_y, 
             head_width=0.12, head_length=0.15,
             fc='orange', ec='darkorange', lw=3, zorder=5)
    
    # N→C arrow in center (keep the arrow)
    ax.annotate('', xy=(0, -0.3), xytext=(0, 0.3),
                arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    # N and C labels at first and last residues (not in center)
    # N at first residue (index 0) - slightly offset from the residue
    n_offset_dist = 0.2
    n_angle = angles[0]
    n_x = (radius + n_offset_dist) * math.cos(n_angle)
    n_y = (radius + n_offset_dist) * math.sin(n_angle)
    ax.text(n_x, n_y, 'N', ha='center', va='center', 
            fontsize=11, color='red', weight='bold')
    
    # C at last residue (index -1) - slightly offset from the residue
    c_offset_dist = 0.22
    c_angle = angles[-1]
    c_x = (radius + c_offset_dist) * math.cos(c_angle)
    c_y = (radius + c_offset_dist) * math.sin(c_angle)
    ax.text(c_x, c_y, 'C', ha='center', va='center', 
            fontsize=11, color='red', weight='bold')
    
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=16, weight='bold', pad=20)
    
    plt.tight_layout()
    return fig

# Main execution
def analyze_amphipathic_helices(seq, window=18, output_dir="/mnt/user-data/outputs", save_files=True, show_plot=False):
    """
    Analyze a protein sequence for amphipathic helices.
    
    Parameters:
    -----------
    seq : str
        Amino acid sequence to analyze
    window : int
        Window size for helix analysis (default: 18)
    output_dir : str
        Directory to save output files (default: "/mnt/user-data/outputs")
    save_files : bool
        Whether to save CSV and plot files (default: True)
    show_plot : bool
        Whether to display the plot (default: False)
    
    Returns:
    --------
    tuple : (best_sequence, max_moment, net_charge)
        best_sequence : str - The sequence with highest hydrophobic moment
        max_moment : float - The maximum hydrophobic moment value
        net_charge : int - The net charge of the best sequence
    """
    global WINDOW
    WINDOW = window
    
    # Analyze sequence
    df = analyze_sequence(seq)
    
    if save_files:
        import os
        os.makedirs(output_dir, exist_ok=True)
        csv_path = os.path.join(output_dir, "heliquest_results.csv")
        df.to_csv(csv_path, index=False)
        print(f"✓ Results saved to '{csv_path}'")
    
    # Find most amphipathic helix
    best_idx = df["HydrophobicMoment"].idxmax()
    best_row = df.loc[best_idx]
    
    # Extract return values
    best_sequence = best_row["Sequence"]
    max_moment = float(best_row["HydrophobicMoment"])
    net_charge = int(best_row["NetCharge"])
    
    if save_files:
        # Print summary
        print("\n" + "=" * 85)
        print("MOST AMPHIPATHIC HELIX")
        print("=" * 85)
        print(f"Helix Number:        #{int(best_row['Helix'])}")
        print(f"Position:            {int(best_row['Start'])}-{int(best_row['End'])}")
        print(f"Sequence:            {best_sequence}")
        print(f"Hydrophobicity (H):  {best_row['Hydrophobicity']:.3f}")
        print(f"Hydrophobic Moment:  {max_moment:.3f}")
        print(f"Net Charge (z):      {net_charge}")
        print(f"Freq Polar:          {best_row['FreqPolar']:.3f}")
        print(f"Freq NonPolar:       {best_row['FreqNonPolar']:.3f}")
        print("=" * 85)
        
        # Create and save wheel plot
        fig = plot_wheel(best_row["Sequence"], best_row["vx"], best_row["vy"],
                         title=f"Most Amphipathic Helix (#{int(best_row['Helix'])})")
        plot_path = os.path.join(output_dir, "heliquest_wheel.png")
        plt.savefig(plot_path, dpi=150, bbox_inches='tight')
        print(f"✓ Wheel plot saved to '{plot_path}'")
        
        if show_plot:
            plt.show()
        else:
            plt.close(fig)
    
    return best_sequence, max_moment, net_charge


if __name__ == "__main__":
    # Example usage
    seq = "SPLTDSGGYYYPKGELDSLLEDIDSPFSIISSQKGKKSEPKRGFPSARFVIDDDSPVAQGKPKKRFLGLF"
    
    # Run analysis
    best_seq, moment, charge = analyze_amphipathic_helices(
        seq, 
        window=18, 
        output_dir="/mnt/user-data/outputs",
        save_files=True,
        show_plot=False
    )
    
    print(f"\nReturned values:")
    print(f"  Best sequence: {best_seq}")
    print(f"  Max moment:    {moment}")
    print(f"  Net charge:    {charge}")
    
    print("=" * 85)