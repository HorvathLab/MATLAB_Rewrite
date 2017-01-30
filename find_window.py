import read_data

method = 'adaptive'
number_of_wdw = 10
chrm_len = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
            159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
            115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
            59128983, 63025520, 48129895, 51304566, 155270560]

def position_window(vaf):
    wdw_variants_num = [(len(vaf) // number_of_wdw)] * (number_of_wdw)
    # wdw_variants_num.append(chrm_variants_num - (number_of_wdw - 1) * (chrm_variants_num//number_of_wdw))
    return wdw_variants_num

def data_window(chrm, vaf):
    start_pos = vaf['POS'].min()
    end_pos = vaf['POS'].max()
    segment_length = round(chrm_len[chrm - 1] / 10)
    wdw_variants_num = list()

    pos = start_pos
    while pos < end_pos:
        segment_vafs = vaf[vaf['POS'] >= pos]
        pos += segment_length
        segment_vafs = segment_vafs[segment_vafs['POS'] < pos]
        wdw_variants_num.append(len(segment_vafs))
    print wdw_variants_num
    return wdw_variants_num

def adaptive_window(chrm, vaf):
    chrm_variants_num = len(vaf) // number_of_wdw
    data_wdw_variants_num = data_window(chrm, vaf)
    wdw_variants_num = list()
    i = 0
    while i < number_of_wdw:
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
    print wdw_variants_num
    return wdw_variants_num

vaf = read_data.getData()
position_window(vaf)
data_window(1, vaf)
adaptive_window(1, vaf)
