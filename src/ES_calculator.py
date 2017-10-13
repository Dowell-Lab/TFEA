__author__ = 'Jonathan Rubin'

import math
import numpy as np
import matplotlib.pyplot as plt

def run(ranked_file,figuredir):
    H = 1500
    ES = list()
    colors = list()
    with open(ranked_file) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            val = float(line[-1])
            if val > H:
                colors.append(0)
            else:
                colors.append(math.fabs(val-H)/H)

    ind = np.arange(0,len(colors))
    ticks = list()
    for val in colors:
        if val > 0:
            ticks.append(1)
        else:
            ticks.append(0)
    F = plt.figure(figsize=(300,20))
    # cbar = plt.colorbar(colors)
    plt.scatter(ind,colors,edgecolor="")
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off')
    plt.savefig(figuredir + 'test.png')

    plt.close()
