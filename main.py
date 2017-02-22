import read_data
import chromosome
import find_window
import KS_test
import scipy


vaf_list = ['nex', 'ntr', 'tex', 'ttr']

def main():
    chrm = find_window.assign_window(17, read_data.getAllVariants(), 'adaptive')

    KS_test.distribution_chromosome(chrm, 'tex')
    KS_test.mode_difference(chrm, 'tex')

    KS_test.distribution_chromosome(chrm, 'ttr')
    KS_test.mode_difference(chrm, 'ttr')

    group_variants_tex = chrm.join_groups('tex')

    texGroup0_tex = [v.tex for v in group_variants_tex[0]]
    texGroup1_tex = [v.tex for v in group_variants_tex[1]]
    chrm.p_values['texGroup_tex'] = scipy.stats.ks_2samp(texGroup0_tex, texGroup1_tex).pvalue

    group_variants_ttr = chrm.join_groups('ttr')

    ttrGroup0_ttr = [v.ttr for v in group_variants_ttr[0]]
    ttrGroup1_ttr = [v.ttr for v in group_variants_ttr[1]]
    chrm.p_values['ttrGroup_ttr'] = scipy.stats.ks_2samp(ttrGroup0_ttr, ttrGroup1_ttr).pvalue

    mode_tex = [w.mode_tex for w in chrm.windows]
    mode_ttr = [w.mode_tex for w in chrm.windows]
    allVariants = chrm.getAllVariants

    if chrm.p_values['texGroup_tex'] > 0.05 or (max(mode_tex) == min(mode_tex)):
        chrm.p_values['texGroup_tex_nex'] = scipy.stats.ks_2samp([v.tex for v in allVariants], [v.nex for v in allVariants]).pvalue
        chrm.p_values['texGroup_tex_ntr'] = scipy.stats.ks_2samp([v.tex for v in allVariants], [v.ntr for v in allVariants]).pvalue
    else:
        texGroup0_nex = [v.nex for v in group_variants_tex[0]]
        chrm.p_values['texGroup0_tex_nex'] = scipy.stats.ks_2samp(texGroup0_tex, texGroup0_nex).pvalue

        texGroup0_ntr = [v.ntr for v in group_variants_tex[0]]
        chrm.p_values['texGroup0_tex_ntr'] = scipy.stats.ks_2samp(texGroup0_tex, texGroup0_ntr).pvalue

        texGroup1_nex = [v.nex for v in group_variants_tex[1]]
        chrm.p_values['texGroup1_tex_nex'] = scipy.stats.ks_2samp(texGroup1_tex, texGroup1_nex).pvalue

        texGroup1_ntr = [v.ntr for v in group_variants_tex[1]]
        chrm.p_values['texGroup1_tex_ntr'] = scipy.stats.ks_2samp(texGroup1_tex, texGroup1_ntr).pvalue

    if chrm.p_values['ttrGroup_ttr'] > 0.05 or (max(mode_ttr) == min(mode_ttr)):
        chrm.p_values['ttrGroup_ttr_nex'] = scipy.stats.ks_2samp([v.ttr for v in allVariants], [v.nex for v in allVariants]).pvalue
        chrm.p_values['ttrGroup_ttr_ntr'] = scipy.stats.ks_2samp([v.ttr for v in allVariants], [v.ntr for v in allVariants]).pvalue
    else:
        ttrGroup0_nex = [v.nex for v in group_variants_ttr[0]]
        chrm.p_values['ttrGroup0_ttr_nex'] = scipy.stats.ks_2samp(ttrGroup0_ttr, ttrGroup0_nex).pvalue

        ttrGroup0_ntr = [v.ntr for v in group_variants_ttr[0]]
        chrm.p_values['ttrGroup0_ttr_ntr'] = scipy.stats.ks_2samp(ttrGroup0_ttr, ttrGroup0_ntr).pvalue

        ttrGroup1_nex = [v.nex for v in group_variants_ttr[1]]
        chrm.p_values['ttrGroup1_ttr_nex'] = scipy.stats.ks_2samp(ttrGroup1_ttr, ttrGroup1_nex).pvalue

        ttrGroup1_ntr = [v.ntr for v in group_variants_ttr[1]]
        chrm.p_values['ttrGroup1_ttr_ntr'] = scipy.stats.ks_2samp(ttrGroup1_ttr, ttrGroup1_ntr).pvalue

    print chrm.p_values


if __name__ == '__main__':
    main()
