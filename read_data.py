import pandas

def getData():
    vaf = pandas.read_csv('001-R.csv')
    return vaf


getData()
