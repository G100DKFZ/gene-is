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
    pos1= values[3]
    M1= cig1.index ("M")
    S1= cig1.index ("S")
    pos1=int (pos1)
    if M1 < S1 :
	spanM1=cig1[:M1]
	spanM1= int(spanM1)
	spanM1=spanM1-1
	correctedMpos=pos1+spanM1
	
        
        print i.rstrip('\n'), correctedMpos
