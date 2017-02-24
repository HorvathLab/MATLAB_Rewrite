import numpy as np
import emd
import matplotlib.pyplot as plt


def plot_distribution_tex(sample):
    np.random.seed(0)

    ideal_3367 = emd.ideal_generator_tex('3367')
    ideal_2575 = emd.ideal_generator_tex('2575')
    ideal_2080 = emd.ideal_generator_tex('2080')
    bins = np.linspace(0, 1, num=31)

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(9, 9))

    ax0.hist(ideal_3367, bins, normed=1, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5)
    ax0.hist(sample, bins, normed=1, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5)

    # Create a histogram by providing the bin edges (unequally spaced).
    ax1.hist(ideal_2575, bins, normed=1, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5)
    ax1.hist(sample, bins, normed=1, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5)

    ax2.hist(ideal_2080, bins, normed=1, histtype='bar', facecolor='r', rwidth=0.8, alpha=0.5)
    ax2.hist(sample, bins, normed=1, histtype='bar', facecolor='b', rwidth=0.8, alpha=0.5)

    fig.tight_layout()
    plt.show()


def plot_grid():
    fig, ax = plt.subplots()
    grid = np.random.rand(4, 4)
    ax.imshow(grid, interpolation=None, cmap='viridis')
    ax.set_title(None)
    plt.show()
    return

# plot_grid()
