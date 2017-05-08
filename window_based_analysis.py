import read_data
import find_window
import pdb
import emd
import numpy
import scipy.stats
from itertools import combinations
import operator

bonf_corr = 4370

class emd_obj:
    def __init__(self, index, emd, fdr):
        self.index = index
        self.emd = emd
        self.fdr = fdr

def get_all_emd(grouped_chrm):
    for comb in combs:
        s1 = [v.tex for v in grouped_chrm.windows[comb[0]].variants]
        s2 = [v.tex for v in grouped_chrm.windows[comb[1]].variants]
        s1 = [abs(v - 0.5) for v in s1]
        s2 = [abs(v - 0.5) for v in s2]
        if s1 != [] and s2 != []:
            p = scipy.stats.ks_2samp(s1, s2).pvalue
            emd_object = emd_obj(comb, emd.emd(emd.cumsum(s1), emd.cumsum(s2)), p * bonf_corr)
        else:
            emd_object = emd_obj(comb, None, None)
        all_emd.append(emd_object)

chrm = find_window.assign_window(17, read_data.getAllVariants(), 'adaptive')
grouped_chrm = chrm
num_window = len(chrm.windows)
combs = list(combinations([x for x in range(num_window)], 2))
all_emd = list()
get_all_emd(chrm)
i = 0
while max([e.fdr for e in all_emd if e.fdr != None]) > 0.05:
    i += 1
    all_emd = [e for e in all_emd if e.emd != None]
    all_emd.sort(key=lambda x: x.emd)
    grouped_chrm.join_windows(all_emd[0].index[0], all_emd[0].index[1])
    all_emd = list()
    get_all_emd(chrm)
# pdb.set_trace()
mins = list()
for win in grouped_chrm.windows:
    bins = numpy.linspace(0.5, 1, num=51)
    if len(win.variants) > 0:
        find_min = list()
        for i in bins:
            ideal = emd.ideal_generator(i)
            emd_val = emd.emd(emd.cumsum([abs(v.tex-0.5) for v in win.variants]), emd.cumsum([abs(v-0.5) for v in ideal]))
            find_min.append(emd_obj(i, emd_val, None))
        print [v.emd for v in find_min]
        find_min.sort(key=lambda x: x.emd)
        mins.append(find_min[0].index)

print mins
