import numpy

bins = 100
max_emd = 1

def ideal_generator_tex(dis_type):
    if dis_type == '0100':
        ideal = numpy.random.beta(0.4, 0.4, 20000)
    else:
        if dis_type == '4555':
            mu_small = 0.45
            mu_large = 0.55
        elif dis_type == '3367':
            mu_small = 0.33
            mu_large = 0.67
        elif dis_type == '2575':
            mu_small = 0.25
            mu_large = 0.75
        elif dis_type == '2080':
            mu_small = 0.2
            mu_large = 0.8
        ideal_small = numpy.random.normal(mu_small, 0.075, 10000)
        ideal_large = numpy.random.normal(mu_large, 0.075, 10000)
        ideal = numpy.append(ideal_small, ideal_large)
        ideal = [k for k in ideal if abs(k - 0.5) < 0.5 ]
    return ideal

def ideal_generator_ttr(dis_type):
    if dis_type == '01b':
        ideal = numpy.random.beta(0.4, 0.4, 20000)
    elif dis_type == '01n':
        ideal_small = numpy.random.normal(0, 0.15, 10000)
        ideal_large = numpy.random.normal(1, 0.15, 10000)
        ideal = numpy.append(ideal_small, ideal_large)
        ideal = [k for k in ideal if abs(k - 0.5) < 0.5 ]
    elif dis_type == '5050':
        ideal = numpy.random.normal(0.5, 0.075, 20000)
        ideal = [k for k in ideal if abs(k - 0.5) < 0.5 ]
    print len(ideal)
    return ideal


def cumsum(data):
    counts = numpy.histogram(data, bins=bins, range=(0, 1))[0]
    length = len(data)
    return numpy.cumsum([float(h)/length for h in counts])

def emd(sample, ideal):
    sub = [abs(x - y) for x, y in zip(sample, ideal)]
    emd = sum(sub) * (1 / float(bins)) / float(max_emd)
    return emd
