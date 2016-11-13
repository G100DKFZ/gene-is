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
    M1= cig1.index ("M")
    S1= cig1.index ("S")
    if S1 < M1 :
        print i
