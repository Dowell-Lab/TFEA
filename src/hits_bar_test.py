import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

def run():
    a = np.random.rand(60000)
    hits = [1 if x>0.5 else 0 for x in a]
    hitlength = len(hits)
    width = hitlength/60000
    # YlOrRd = plt.get_cmap('YlOrRd')
    # hist,_=np.histogram(hits,bins=width)
    # newhits = [0]*len(hits)
    # for i in range(0,len(hits),window):
    #     if sum(hits[i:i+window]) > 1:
    #         newhits[i+window/2] = 1

    newhits = list()
    for i in range(0,hitlength,width):
        newhits.append(float(sum(hits[i:i+width])))

    newhitmax = float(max(newhits))
    newhitlength = len(newhits)

    alphas = [(x/newhitmax) for x in newhits]
    rgba_colors = np.zeros((newhitlength,4))
    rgba_colors[:,0] = 0.0
    rgba_colors[:,1] = 0.0
    rgba_colors[:,2] = 0.0
    rgba_colors[:,3] = alphas

    F = plt.figure(figsize=(15,6))
    ax1 = plt.subplot(111)
    print range(1,hitlength+1,width)[-10:], ([1]*newhitlength)[:10], rgba_colors[:10]
    ax1.bar(range(1,hitlength+1,width), [1]*newhitlength , color=rgba_colors,edgecolor = "none",width=width)
    ax1.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left='off',        # ticks along the bottom edge are off
        right='off',       # ticks along the top edge are off
        labelleft='off')
    ax1.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on')
    ax1.set_xlim([0,hitlength])
    ax1.set_ylim([0,1])
    ax1.set_xlabel('Rank in Ordered Dataset', fontsize=14)
    ax1.set_ylabel('Hits', fontsize=14)
    # plt.savefig(figuredir + MOTIF_FILE + '_enrichment_plot.svg')
    plt.show()

if __name__ == "__main__":
    run()