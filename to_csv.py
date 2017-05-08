import csv

e = list()
f = open('001-R.txt', 'r')
for line in f:
    all_elements = line[:-1].split('\t')
    if all_elements[0] == '17':
        e.append(all_elements)
with open('001-R.csv', "w") as output:
    writer = csv.writer(output, dialect='excel')
    writer.writerows(e)
