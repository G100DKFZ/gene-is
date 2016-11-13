import csv
import sys

Lines= sys.argv[1]
datain=sys.argv[2]


Lines= open(Lines, "r")
for x in Lines:
    x=x
n=int(x)
a=[""]*n
j=0
datain=open(datain, "r")
for i in datain:
    values=i.split()
    cig1= values[5]
    if cig1.count("M")==1 and "S" not in cig1:
        print i
