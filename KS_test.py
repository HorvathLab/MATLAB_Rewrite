import find_window
import read_data
import numpy as np
from collections import Counter
import scipy.stats

def get_distribution(vaf_wdw, vaf):
    variant_in_wdw = vaf_wdw.variants
    vaf_distance = [round(abs(getattr(v, vaf)  - 0.5), 2) for v in variant_in_wdw]
    bins = np.linspace(0, 0.5, num = 51)
    hist = np.histogram(vaf_distance, bins=bins)
    idx = np.argmax(hist[0])
    return round(hist[1][idx], 2)

def distribution_chromosome(chrm, vaf):
    all_windows = chrm.windows
    for w in all_windows:
        if vaf == 'ttr' and get_distribution(w, vaf) == 0.34:
            setattr(w, 'mode_' + vaf, 0)
        else:
            setattr(w, 'mode_' + vaf, get_distribution(w, vaf))
    return

def mode_difference(chrm, vaf):
    max_dis = [getattr(w, 'mode_' + vaf) for w in chrm.windows]
    maxr = max(max_dis)
    minr = min(max_dis)
    for w in chrm.windows:
        if getattr(w, 'mode_' + vaf) <= (maxr + minr)/2:
            w.group = 0
            setattr(w, 'group_' + vaf, 0)
        else:
            setattr(w, 'group_' + vaf, 1)
    return
