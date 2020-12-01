# This script is adapted from a tutorial by the James Hutton Institute, original creator is unknown.
# URL: https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/05-blast_for_rbh.html

"""
RBH_plotting (RBH_tools).
Version: 1.0
Last modified: 16/11/2020
Author: Cailean Carter
"""

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
import seaborn as sns
import numpy as np
import sys


def plot(filepath):

    """Plots two graphs:
    - A heatmap of query and subject coverage.
    - A distribution plot of RBH normalised bit-scores.

    Takes the output tabular file from BLAST_RBH tool on Galaxy.
    Takes file path to the tabular file as positional argument.
    """

    rbh = pd.read_csv(filepath, sep="\t")

    #normalise bitscore
    rbh["norm_bitscore"] = rbh.bitscore/rbh.length

    # Plot 2D density histograms

    # Calculate 2D density histograms for counts of matches at several coverage levels
    (H, xedges, yedges) = np.histogram2d(rbh.A_qcovhsp, rbh.B_qcovhsp, bins=20)

    # Create a 1x1 figure array
    fig, ax = plt.subplots(1, 2, figsize=(6, 6))

    # Plot histogram for RBBH
    his = ax[0].imshow(H, cmap=plt.cm.Blues, norm=LogNorm(),
                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                    origin='lower', aspect=1)
    ax[0].set_title("RBH")
    ax[0].set_xlabel("Query")
    ax[0].set_ylabel("Subject")

    # Add colourbar
    fig.colorbar(his, ax=ax[0])

    # Plot distribution of RBH bitscores
    sns.distplot(rbh.norm_bitscore, color="b", axlabel="RBH normalised bitscores", ax=ax[1])

    plt.show()

if __name__ == "__main__":
    plot(sys.argv[1])