import numpy

max_emd = 0.5
num_bins = 101
bins = numpy.linspace(0, 0.5, num=num_bins)
bin_width = bins[2] - bins[1]

def ideal_generator(dis_type):
    mu_small = 1 - dis_type
    mu_large = dis_type
    ideal_small = numpy.random.normal(mu_small, 0.075, 10000)
    ideal_large = numpy.random.normal(mu_large, 0.075, 10000)
    ideal = numpy.append(ideal_small, ideal_large)
    ideal = [k for k in ideal if abs(k - 0.5) < 0.5 ]
    return ideal

def cumsum(data):
    if data:
        counts = numpy.histogram(data, bins=bins)[0]
        length = len(data)
        return numpy.cumsum([float(h)/length for h in counts])
    # else:
    #     raise RuntimeError

def emd(sample1, sample2):
    sub = [abs(x - y) for x, y in zip(sample1, sample2)]
    emd = sum(sub) * bin_width / float(max_emd)
    return emd
