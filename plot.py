import numpy as np
import emd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def plot_distribution_tex(sample):
    np.random.seed(0)

    ideal_3367 = emd.ideal_generator_tex('3367')
    ideal_2575 = emd.ideal_generator_tex('2575')
    ideal_2080 = emd.ideal_generator_tex('2080')
    bins = np.linspace(0, 1, num=31)
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(9, 9))

    ax0.hist(ideal_3367, bins, normed=0, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5, weights=np.ones_like(ideal_3367)/float(len(ideal_3367)))
    ax0.hist(sample, bins, normed=0, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5, weights=np.ones_like(sample)/float(len(sample)))

    ax1.hist(ideal_2575, bins, normed=0, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5, weights=np.ones_like(ideal_2575)/float(len(ideal_2575)))
    ax1.hist(sample, bins, normed=0, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5, weights=np.ones_like(sample)/float(len(sample)))

    ax2.hist(ideal_2080, bins, normed=0, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5, weights=np.ones_like(ideal_2080)/float(len(ideal_2080)))
    ax2.hist(sample, bins, normed=0, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5, weights=np.ones_like(sample)/float(len(sample)))

    fig.tight_layout()
    plt.savefig('samp_1_chrm_17_dis.png')


def plot_grid(x, y):
    plt.hist2d(x, y, bins=np.linspace(0, 1, num=31))
    plt.show()

    return

# plot_grid()
