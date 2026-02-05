import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



# ==============================
# PARAMETERS
# ==============================
WINDOW_SIZE = 18
ANGLE_DEG = 100  # rotation per residue in alpha helix

# ==============================
# Eisenberg Hydrophobicity Scale
# ==============================
HYDRO = {
    'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90,
    'C': 0.29, 'Q': -0.85, 'E': -0.74, 'G': 0.48,
    'H': -0.40, 'I': 1.38, 'L': 1.06, 'K': -1.50,
    'M': 0.64, 'F': 1.19, 'P': 0.12, 'S': -0.18,
    'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
}

POLAR = set("STNQHYG")
NONPOLAR = set("AVLIMFWCP")


# ==============================
# CORE CALCULATIONS
# ==============================

def mean_hydrophobicity(segment):
    return sum(HYDRO.get(aa, 0) for aa in segment) / len(segment)


def hydrophobic_moment(segment):
    angle = math.radians(ANGLE_DEG)
    x, y = 0.0, 0.0

    for i, aa in enumerate(segment):
        h = HYDRO.get(aa, 0)
        theta = i * angle
        x += h * math.cos(theta)
        y += h * math.sin(theta)

    return math.sqrt(x*x + y*y) / len(segment)


def net_charge(segment):
    charge = 0
    for aa in segment:
        if aa in "KR":
            charge += 1
        elif aa in "DE":
            charge -= 1
    return charge


def residue_frequencies(segment):
    polar = sum(aa in POLAR for aa in segment)
    nonpolar = sum(aa in NONPOLAR for aa in segment)
    n = len(segment)
    return polar / n, nonpolar / n



# ==============================
# Residue Color Mapping
# ==============================
def residue_color(aa):
    if aa in "KR":
        return "blue"
    elif aa in "DE":
        return "red"
    elif aa in "STNQ":
        return "purple"
    elif aa == "G":
        return "lightgrey"
    elif aa == "P":
        return "limegreen"
    else:
        return "yellow"


# ==============================
# HELIQUEST-STYLE HELICAL WHEEL
# ==============================
def plot_heliquest_wheel(segment, title="Helical Wheel"):
    n = len(segment)
    angle_step = np.deg2rad(100)
    radius = 1.2

    angles = np.arange(n) * angle_step
    x = radius * np.cos(angles)
    y = radius * np.sin(angles)

    fig, ax = plt.subplots(figsize=(7,7))

    # Draw backbone circle
    circle = plt.Circle((0,0), radius, color='black', fill=False, linewidth=2)
    ax.add_patch(circle)

    # Draw residue connections (helix trace)
    for i in range(n):
        j = (i + 3) % n  # connect i to i+3 (helical turn)
        ax.plot([x[i], x[j]], [y[i], y[j]], color="black", linewidth=1)

    # Draw residues as colored circles
    for i, aa in enumerate(segment):
        ax.scatter(x[i], y[i], s=1600, color=residue_color(aa), edgecolors='black', zorder=3)
        ax.text(x[i], y[i], aa, ha='center', va='center',
                fontsize=14, weight='bold',
                color='white' if residue_color(aa) in ["blue","red","purple"] else 'black')

    # Direction arrow (N→C)
    # ==============================
    # Direction Arrow (HeliQuest Style)
    # ==============================
    arrow_len = radius * 0.9

    # Vertical arrow pointing down
    ax.arrow(0, arrow_len/2, 0, -arrow_len,
            head_width=0.15, head_length=0.2,
            fc='black', ec='black', linewidth=2)

    # N label at top of arrow
    ax.text(0, arrow_len/2 + 0.2, "N",
            ha='center', va='bottom', fontsize=14, color='red', weight='bold')


    ax.set_title(title, fontsize=16)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.show()


# ==============================
# MAIN ANALYSIS FUNCTION
# ==============================

def heliquest_analysis(sequence, window=WINDOW_SIZE, make_plot=False):
    sequence = sequence.upper()
    rows = []

    for i in range(len(sequence) - window + 1):
        segment = sequence[i:i+window]

        H = mean_hydrophobicity(segment)
        muH = hydrophobic_moment(segment)
        z = net_charge(segment)
        f_polar, f_nonpolar = residue_frequencies(segment)

        rows.append({
            "Helix": i + 1,
            "Start": i + 1,
            "End": i + window,
            "Sequence": segment,
            "Hydrophobicity": round(H, 3),
            "HydrophobicMoment": round(muH, 3),
            "NetCharge": z,
            "FreqPolar": round(f_polar, 3),
            "FreqNonPolar": round(f_nonpolar, 3)
        })

    df = pd.DataFrame(rows)

    if make_plot:
        best = df["HydrophobicMoment"].idxmax()
        plot_heliquest_wheel(df.loc[best, "Sequence"],
                           title=f"Most Amphipathic Helix (#{df.loc[best,'Helix']})")

    return df


# ==============================
# RUN AS SCRIPT
# ==============================

if __name__ == "__main__":
    sequence = "SPLTDSGGYYYPKGELDSLLEDIDSPFSIISSQKGKKSEPKRGFPSARFVIDDDSPVAQGKPKKRFLGLF"

    df = heliquest_analysis(sequence, make_plot=True)
    print(df.head())

    df.to_csv("heliquest_pandas_results.csv", index=False)
    print("\nSaved to heliquest_pandas_results.csv")
