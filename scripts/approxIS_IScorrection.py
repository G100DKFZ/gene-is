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
    pos1=values[3]
    cig1=int(cig1)
    pos1=int(pos1)
    corr=cig1+pos1
    corr1=corr-1
    print corr1
