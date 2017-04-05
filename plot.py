import numpy as np
import emd
import find_window
import read_data
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def plot_distribution_tex(sample):
    np.random.seed(0)

    ideal_3367 = emd.ideal_generator_tex('3367')
    ideal_2575 = emd.ideal_generator_tex('2575')
    ideal_2080 = emd.ideal_generator_tex('2080')
    bins = np.linspace(0, 1, num=31)
    fig = plt.figure()
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(9, 9))

    ax0.hist(ideal_3367, bins, normed=0, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5, weights=np.ones_like(ideal_3367)/float(len(ideal_3367)))
    ax0.hist(sample, bins, normed=0, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5, weights=np.ones_like(sample)/float(len(sample)))

    ax1.hist(ideal_2575, bins, normed=0, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5, weights=np.ones_like(ideal_2575)/float(len(ideal_2575)))
    ax1.hist(sample, bins, normed=0, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5, weights=np.ones_like(sample)/float(len(sample)))

    ax2.hist(ideal_2080, bins, normed=0, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5, weights=np.ones_like(ideal_2080)/float(len(ideal_2080)))
    ax2.hist(sample, bins, normed=0, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5, weights=np.ones_like(sample)/float(len(sample)))

    fig.tight_layout()
    plt.savefig('samp_1_chrm_17_dis.png')
    plt.figure()
    return


def plot_grid(x, y):
    fig2 = plt.hist2d(x, y, bins=np.linspace(0, 1, num=31))
    plt.savefig('tex_ttr_bimodal.png')
    plt.figure()
    return


def plot_chrm(chrm):
    allVariants = chrm.getAllVariants()
    pos = [v.pos for v in allVariants]
    tex = [v.tex for v in allVariants]
    ttr = [v.ttr for v in allVariants]
    fig = plt.figure()
    plt.scatter(pos, tex, c='red', s=0.5, alpha=1)
    plt.scatter(pos, ttr, c='orange', s=0.5, alpha=1)
    plt.show()


# chrm = find_window.assign_window(17, read_data.getAllVariants(), 'adaptive')
# plot_chrm(chrm)
