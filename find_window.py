import chromosome

method = 'adaptive'
number_of_wdw = 10
chrm_len = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
            159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
            115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
            59128983, 63025520, 48129895, 51304566, 155270560]

def position_window(chrm_variants_num):
    wdw_variants_num = [(chrm_variants_num // number_of_wdw)] * (number_of_wdw - 1)
    wdw_variants_num.append(chrm_variants_num - (number_of_wdw - 1) * (chrm_variants_num//number_of_wdw))
    return wdw_variants_num

def data_window(chrom, allVariants):
    start_pos = min([v.pos for v in allVariants])
    end_pos = max([v.pos for v in allVariants])
    segment_length = round(chrm_len[chrom - 1] / number_of_wdw)
    wdw_variants_num = list()

    pos = start_pos
    while pos < end_pos:
        segment_vafs = [v for v in allVariants if v.pos >= pos]
        pos += segment_length
        segment_vafs = [v for v in segment_vafs if v.pos < pos]
        wdw_variants_num.append(len(segment_vafs))
    return wdw_variants_num

def adaptive_window(chrom, allVariants):
    chrm_variants_num = len(allVariants) // number_of_wdw
    data_wdw_variants_num = data_window(chrom, allVariants)
    wdw_variants_num = list()
    i = 0
    while i < number_of_wdw:
        if wdw_variants_num and wdw_variants_num[-1] < chrm_variants_num:
            wdw_variants_num[-1] += data_wdw_variants_num[i]
            i += 1
        else:
            if data_wdw_variants_num[i] < chrm_variants_num:
                if i != number_of_wdw - 1:
                    wdw_variants_num.append(data_wdw_variants_num[i] + data_wdw_variants_num[i + 1])
                    i += 2
                elif i == number_of_wdw - 1:
                    wdw_variants_num[-1] += data_wdw_variants_num[i]
                    i += 1
            else:
                wdw_variants_num.append(data_wdw_variants_num[i])
                i += 1
    return wdw_variants_num

def assign_window(chrom, allVariants, type):
    chrm = chromosome.chromosome()
    if type == 'data':
        wdw_variants_num = data_window(chrom, allVariants)
    elif type == 'position':
        wdw_variants_num = position_window(len(allVariants))
    elif type == 'adaptive':
        wdw_variants_num = adaptive_window(chrom, allVariants)
    else:
        raise ValueError('%s is not invalid. Please enter data/position/adaptive as type' % type)
    start = 0
    for w in wdw_variants_num:
        window = chromosome.window()
        end = start + w
        window.variants = allVariants[start:end]
        chrm.windows.append(window)
        start += w
    return chrm
