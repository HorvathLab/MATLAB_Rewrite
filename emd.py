import numpy

bins = 100
max_emd = 1

def ideal_generator(dis_type):
    mu_small = 1 - dis_type
    mu_large = dis_type
    ideal_small = numpy.random.normal(mu_small, 0.075, 10000)
    ideal_large = numpy.random.normal(mu_large, 0.075, 10000)
    ideal = numpy.append(ideal_small, ideal_large)
    ideal = [k for k in ideal if abs(k - 0.5) < 0.5 ]
    return ideal

def cumsum(data):
    counts = numpy.histogram(data, bins=bins, range=(0, 1))[0]
    length = len(data)
    return numpy.cumsum([float(h)/length for h in counts])

def emd(sample, ideal):
    sub = [abs(x - y) for x, y in zip(sample, ideal)]
    emd = sum(sub) * (1 / float(bins)) / float(max_emd)
    return emd
