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
    seq1= values[9]
    M1= cig1.index ("M")
    S1= cig1.index ("S")
    if S1 < M1 :
	span=cig1[:S1]
	span= int(span)
	seq=seq1[:span]
        
        print i.rstrip('\n'), span, seq.rstrip('\n')
