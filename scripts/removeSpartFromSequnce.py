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
    ID= values[0]
    seqFull= values[1]
    seq=values[2]
    
    if seq in seqFull:
	print ID
	print(seqFull.replace(seq,""))

