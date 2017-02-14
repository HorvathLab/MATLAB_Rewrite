import find_window
import read_data
from collections import Counter
import scipy.stats

vaf_list = ['tex', 'ttr', 'nex', 'ntr']

def get_distribution(vaf_wdw, vaf):
    variant_in_wdw = vaf_wdw.variants
    vaf_distance = [round(abs(getattr(v, vaf) - 0.5), 2) for v in variant_in_wdw]
    mode = Counter(vaf_distance).most_common()
    i = 0
    maxima = [mode[0][0]]
    while True:
        if mode[i][1] == mode[i + 1][1]:
            maxima.append(mode[i + 1][0])
            i += 1
        else:
            return min(maxima)
            break

def distribution_chromosome(chrm, vaf):
    all_windows = chrm.windows
    for w in all_windows:
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

allVariants = read_data.getAllVariants()
chrm = find_window.assign_window(allVariants, 'adaptive')
for vaf in vaf_list:
    distribution_chromosome(chrm, vaf)
    mode_difference(chrm, vaf)

group_variants = chrm.join_groups('tex')

data1 = [v.nex for v in group_variants[0]]
data2 = [v.ntr for v in group_variants[0]]
p_value = scipy.stats.ks_2samp(data1, data2)
print p_value
