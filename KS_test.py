import pandas
import read_data
import find_window
import scipy

def get_distribution(vaf_wdw, type):
    type_vaf = pandas.DataFrame(columns=[type])
    type_vaf.loc[:,type] = vaf_wdw.loc[:,type].astype(float)
    type_vaf.loc[:,'folded'] = type_vaf.loc[:,(type)] - 0.5
    type_vaf.loc[:,'folded'] = type_vaf.loc[:,'folded'].abs().round(2)
    # type_vaf = type_vaf.sort_values(by = 'folded', ascending = 1)
    # print type_vaf
    count_vaf = type_vaf.groupby('folded').count()
    return count_vaf[type].idxmax()

vaf = read_data.getData()
wdw_variants_num = find_window.adaptive_window(1, vaf)
# get_distribution(vaf, 'Tex')

def distribution_chromosome(vaf, wdw_variants_num, type):
    start = 0
    end = 0
    max_dis = list()
    for l in wdw_variants_num:
        end = start + l
        vaf_wdw = vaf[start:end]
        start += l
        max_dis.append(get_distribution(vaf_wdw, type))
    return max_dis

Tex_max = distribution_chromosome(vaf, wdw_variants_num, 'Tex')
Ttr_max = distribution_chromosome(vaf, wdw_variants_num, 'Ttr')

def mode_difference(max_dis):
    grouping = [[], []]
    maxr = max(max_dis)
    minr = min(max_dis)
    for index, value in enumerate(max_dis):
        if value <= (maxr + minr)/2:
            grouping[0].append(index+1)
        else:
            grouping[1].append(index+1)
    print grouping
    return grouping

mode_difference(Tex_max)
mode_difference(Ttr_max)
