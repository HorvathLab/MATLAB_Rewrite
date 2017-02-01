import pandas

def getData():
    vaf = pandas.read_csv('001-R.csv')
    vaf = vaf.sort_values(by = 'POS', ascending = 1)
    return vaf


getData()
