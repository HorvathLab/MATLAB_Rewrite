import csv
import chrm


def getAllVariants():
    allVariants = list()
    with open('001-R.csv', 'rb') as csvfile:
        allData = csv.reader(csvfile)
        allData.next()
        for line in allData:
            variant = chrm.variant(line[0], int(line[1]), line[2], line[3], float(line[4]),
                                   float(line[5]), float(line[6]), float(line[7]))
            allVariants.append(variant)
    return allVariants
