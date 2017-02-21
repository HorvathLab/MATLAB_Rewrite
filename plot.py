import numpy as np
import emd
import matplotlib.pyplot as plt


def plot_distribution():
    np.random.seed(0)

    ideal = emd.ideal_generator('3367')
    bins = np.linspace(0, 1, num=31)

    fig, (ax0, ax1) = plt.subplots(nrows=3, figsize=(9, 9))

    ax0.hist(ideal, 0.2, normed=1, histtype='bar', facecolor='g', rwidth=0.8)
    ax0.set_title('stepfilled')

    # Create a histogram by providing the bin edges (unequally spaced).
    ax1.hist(ideal, bins, normed=1, histtype='bar', rwidth=0.8)
    ax1.set_title('unequal bins')

    plt.show()


def plot_grid():
    fig, ax = plt.subplots()
    grid = np.random.rand(4, 4)
    ax.imshow(grid, interpolation=None, cmap='viridis')
    ax.set_title(None)
    plt.show()
    return

plot_grid()

# methods = [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
#            'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
#            'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']
#
# np.random.seed(0)
# grid = np.random.rand(4, 4)
#
# fig, axes = plt.subplots(3, 6, figsize=(12, 6),
#                          subplot_kw={'xticks': [], 'yticks': []})
#
# fig.subplots_adjust(hspace=0.3, wspace=0.05)
#
# for ax, interp_method in zip(axes.flat, methods):
#     ax.imshow(grid, interpolation=interp_method, cmap='viridis')
#     ax.set_title(interp_method)
#
# plt.show()
