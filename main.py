import read_data
import chromosome
import find_window
import KS_test
import scipy
import emd
import plot
import pdb

bonf_corr = 72*22*2
sig_level = 0.05/bonf_corr
vaf_list = ['nex', 'ntr', 'tex', 'ttr']
EMD_ideal = [(51 + x) * 0.01 for x in range(50)]

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

    if 'texGroup_tex_nex' in chrm.p_values.keys():
        if chrm.p_values['texGroup_tex_nex'] < sig_level:
            chrm.bimodal_tex = [x for x in range(len(chrm.windows))]
            chrm.unimodal_tex = []
        else:
            chrm.unimodal_tex = [x for x in range(len(chrm.windows))]
            chrm.bimodal_tex = []
    elif 'texGroup0_tex_nex' in chrm.p_values.keys() and 'texGroup1_tex_nex' in chrm.p_values.keys():
        grouping_tex = [w.group_tex for w in chrm.windows]
        chrm.bimodal_tex = list()
        chrm.unimodal_tex = list()
        if chrm.p_values['texGroup0_tex_nex'] < sig_level:
            chrm.bimodal_tex += [i for i,x in enumerate(grouping_tex) if x == 0]
            chrm.bimodal['tex'].append(0)
        else:
            chrm.unimodal_tex += [i for i,x in enumerate(grouping_tex) if x == 0]

        if chrm.p_values['texGroup1_tex_nex'] < sig_level:
            chrm.bimodal_tex += [i for i,x in enumerate(grouping_tex) if x == 1]
            chrm.bimodal['tex'].append(1)
        else:
            chrm.unimodal_tex += [i for i,x in enumerate(grouping_tex) if x == 1]

        grouping_ttr = [w.group_ttr for w in chrm.windows]
        chrm.bimodal_ttr = list()
        chrm.unimodal_ttr = list()
        if chrm.p_values['ttrGroup0_ttr_ntr'] < sig_level:
            chrm.bimodal_ttr += [i for i,x in enumerate(grouping_ttr) if x == 0]
            chrm.bimodal['ttr'].append(0)
        else:
            chrm.unimodal_ttr += [i for i,x in enumerate(grouping_ttr) if x == 0]

        if chrm.p_values['ttrGroup1_ttr_ntr'] < sig_level:
            chrm.bimodal_ttr += [i for i,x in enumerate(grouping_ttr) if x == 1]
            chrm.bimodal['ttr'].append(1)
        else:
            chrm.unimodal_ttr += [i for i,x in enumerate(grouping_ttr) if x == 1]

        if chrm.bimodal_tex != [] and chrm.bimodal_ttr != []:
            common = list(set(chrm.bimodal_tex).intersection(set(chrm.bimodal_ttr)))
            common_variants = list()
            for w in common:
                common_variants += chrm.windows[w].variants
            chrm.p_values['common_tex_ttr'] = scipy.stats.ks_2samp([v.tex for v in common_variants], [v.ttr for v in common_variants]).pvalue

        if chrm.bimodal_tex != []:
            bi_variants_tex = list()
            for w in chrm.bimodal_tex:
                bi_variants_tex += chrm.windows[w].variants
            chrm.p_values['texGroup_tex_ttr_bimodal'] = scipy.stats.ks_2samp([v.tex for v in bi_variants_tex], [v.ttr for v in bi_variants_tex]).pvalue

        if chrm.bimodal_ttr != []:
            bi_variants_ttr = list()
            for w in chrm.bimodal_ttr:
                bi_variants_ttr += chrm.windows[w].variants
            chrm.p_values['ttrGroup_tex_ttr_bimodal'] = scipy.stats.ks_2samp([v.tex for v in bi_variants_ttr], [v.ttr for v in bi_variants_ttr]).pvalue

        #EMD
        EMD = list()
        for ideal_type in EMD_ideal:
            ideal = emd.ideal_generator(ideal_type)
            EMD.append(emd.emd(emd.cumsum(texGroup0_tex), emd.cumsum(ideal)))
        chrm.emd_group0tex = EMD.index(min(EMD)) + 51

        EMD = list()
        for ideal_type in EMD_ideal:
            ideal = emd.ideal_generator(ideal_type)
            EMD.append(emd.emd(emd.cumsum(texGroup1_tex), emd.cumsum(ideal)))
        chrm.emd_group1tex = EMD.index(min(EMD)) + 51

        EMD = list()
        for ideal_type in EMD_ideal:
            ideal = emd.ideal_generator(ideal_type)
            EMD.append(emd.emd(emd.cumsum(ttrGroup0_ttr), emd.cumsum(ideal)))
            chrm.emd_group0ttr = EMD.index(min(EMD)) + 51

        EMD = list()
        for ideal_type in EMD_ideal:
            ideal = emd.ideal_generator(ideal_type)
            EMD.append(emd.emd(emd.cumsum(ttrGroup1_ttr), emd.cumsum(ideal)))
            chrm.emd_group1ttr = EMD.index(min(EMD)) + 51

        result = open('result.txt', 'w')
        features = ['p_values', 'emd_group0tex', 'emd_group0ttr', 'emd_group1tex', 'emd_group1ttr']
        print chrm.emd_group0tex
        print chrm.emd_group1tex
        print chrm.emd_group0ttr
        print chrm.emd_group1ttr

        # for f in features:
        #     result.write(f + '\n')
        #     chrm_feature = getattr(chrm, f)
        #     for k in chrm_feature.keys():
        #         result.write(str(k) + ': %f\n' % chrm_feature[k])
        # result.close()

        # plot.plot_distribution_tex(texGroup1_tex)
        # texGroup1_ttr = [v.ttr for v in group_variants_tex[1]]
        # plot.plot_grid(texGroup1_ttr, texGroup1_tex)


if __name__ == '__main__':
    main()
