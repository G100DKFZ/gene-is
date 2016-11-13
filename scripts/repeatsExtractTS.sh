vectorString=$1
workingDir=$2
scriptDir=$3
alScore=$4
aligner=$5
refGenomeIndex=$6
minIden=$7
range=$8

 
awk  -v vectorString=$vectorString  '($2!=vectorString)' $workingDir/filtered.bst > $workingDir/tocheckrepeats1.txt

sed 's/@/\t/g'   $workingDir/tocheckrepeats1.txt | awk  -v vectorString=$vectorString  '($2==vectorString)'  | cut -f1  | sort | uniq  > $workingDir/tocheckrepeats1.txt5
grep -Fwf  $workingDir/tocheckrepeats1.txt5 $workingDir/tocheckrepeats1.txt > $workingDir/tocheckrepeats1.txt6

cut -f1 $workingDir/tocheckrepeats1.txt6 | sort | uniq -d | cut -d "@" -f1 | sort | uniq  > $workingDir/tocheckrepeats1.txt7


grep -Fwf $workingDir/tocheckrepeats1.txt7 $workingDir/resultsNoDupSingle.csv | cut -d " " -f1 | sort | uniq  > $workingDir/tocheckrepeats1.txt8
cat $workingDir/resultsNoDupSingle.csv > $workingDir/selected
grep -v -Fwf  $workingDir/tocheckrepeats1.txt8  $workingDir/tocheckrepeats1.txt7 >  $workingDir/tocheckrepeats1.txt9

grep -Fwf $workingDir/tocheckrepeats1.txt9 $workingDir/filtered.bst > $workingDir/tocheckrepeats1.txt10

sort -k1,1 -u $workingDir/tocheckrepeats1.txt10 | awk -v vectorString=$vectorString  '($2==vectorString)' | cut -d "@" -f1 > $workingDir/tocheckrepeats1.txt11

grep -v -Fwf   $workingDir/tocheckrepeats1.txt11  $workingDir/tocheckrepeats1.txt9 >  $workingDir/tocheckrepeats1.txt12
grep -Fwf  $workingDir/tocheckrepeats1.txt12  $workingDir/filtered.bst >  $workingDir/tocheckrepeats1.txt13


sort -k1,1 -k12,12nr $workingDir/tocheckrepeats1.txt13  >  $workingDir/tocheckrepeats1.txt13.z1
sed 's/\t/DDDD/g' $workingDir/tocheckrepeats1.txt13.z1 > $workingDir/tocheckrepeats1.txt13.z2
sort -k1,1  -u $workingDir/tocheckrepeats1.txt13.z1 | sed 's/\t/DDDD/g' > $workingDir/tocheckrepeats1.txt13.z3  
grep -A1 -Fwf  $workingDir/tocheckrepeats1.txt13.z3  $workingDir/tocheckrepeats1.txt13.z2 | sed 's/DDDD/\t/g' |  sed 's/--//g' | sed '/^$/d' > $workingDir/tocheckrepeats1.txt13.z4

sort -k1,1 $workingDir/tocheckrepeats1.txt13.z4 | cut -f1 | uniq -d > $workingDir/tocheckrepeats1.txt13.z5
grep -Fwf  $workingDir/tocheckrepeats1.txt13.z5  $workingDir/tocheckrepeats1.txt13.z4 >  $workingDir/tocheckrepeats1.txt13.z6

awk '{getline b;printf("%s %s\n",$0,b)}'  $workingDir/tocheckrepeats1.txt13.z6 | sed 's/ /\t/g'| cut -f 1,2,3,4,12,14,15,16,24 | awk -v vectorString=$vectorString  '($6!=vectorString)'| cut -f1,2,3,4,5,7,8,9 > $workingDir/tocheckrepeats1.txt13.z7

awk '{ print  $0, ($6/$3),  ($7/$4), ($8/$5)}' $workingDir/tocheckrepeats1.txt13.z7 |  sed 's/ /\t/g' | awk  -v alScore=$alScore '($9>=alScore && $10>=alScore && $11>=alScore )' > $workingDir/tocheckrepeats1.txt13.z8
 cut -f1 $workingDir/tocheckrepeats1.txt13.z8 > $workingDir/tocheckrepeats1.txt13.z9
grep -Fwf $workingDir/tocheckrepeats1.txt13.z9 $workingDir/tocheckrepeats1.txt13.z6 | sort -k1,1 -u > $workingDir/tocheckrepeats1.txt13.z10
mv $workingDir/tocheckrepeats1.txt13.z10   $workingDir/tocheckrepeats1.txt15.filtered.bst 

 sort -k1,1 -k12,12nr $workingDir/tocheckrepeats1.txt15.filtered.bst > $workingDir/tocheckrepeats1.txt15.filtered.sort.bst
awk -F"\t" -v vectorString=$vectorString -f $scriptDir/extractIS.awk $workingDir/tocheckrepeats1.txt15.filtered.sort.bst  | grep $vectorString > $workingDir/repeats.resultsNoDup.csv;
sort -k1,1 -u $workingDir/repeats.resultsNoDup.csv > $workingDir/repeats.resultsNoDupSingle.csv



######################################################################################
######################################################################################
cut -f1 $workingDir/filtered.bst | sed 's/@/\t/g' | sort -k1,1 -u | awk '{ print  ">"$1,$6}' | sed 's/ /\n/g' > $workingDir/filtered.temp.fa
cut -d ' ' -f1 $workingDir/selected | sort -k1,1 -u > $workingDir/selected1
 grep -A1 -Fwf $workingDir/selected1 $workingDir/filtered.temp.fa | sed 's/--//g' | sed '/^\s*$/d' > $workingDir/filtered.temp1.fa
$aligner $refGenomeIndex $workingDir/filtered.temp1.fa -out=blast8 -minIdentity=$minIden $workingDir/filtered.temp.bst
#/home/saira2/bin/x86_64/blat /home/saira2/INDEXES/hg38.v2_vectorSeq.fa.2bit  $workingDir/filtered.temp.fa -out=blast8 -minIdentity=95 $workingDir/filtered.temp.bst

awk  -v vectorString=$vectorString  '($2!=vectorString)' $workingDir/filtered.temp.bst > $workingDir/tocheckrepeats1.temp.txt
sort -k1,1 -k12,12nr $workingDir/tocheckrepeats1.temp.txt > $workingDir/filtered.sort.temp.bst

awk '{n=$1}l1!=n{if(p)print l0; print; p=0}l1==n{p=1}{l0=$0; l1=n}END{print}'  $workingDir/filtered.sort.temp.bst | sort -k1,1 -k12,12n | awk '{if (x[$1]) { x_count[$1]++; print $0; if (x_count[$1] == 1) { print x[$1] } } x[$1] =$0}'  |  awk  '!_[$1]++' >   $workingDir/tocheckrepeats2.temp.txt

awk '{if (x[$1]) { x_count[$1]++; print $0; if (x_count[$1] == 1) { print x[$1] } } x[$1] = $0}' $workingDir/filtered.sort.temp.bst  |  awk  '!_[$1]++'   > $workingDir/tocheckrepeats3.temp.txt

paste  $workingDir/tocheckrepeats2.temp.txt  $workingDir/tocheckrepeats3.temp.txt | cut -f 1,2,3,4,12,15,16,24 |  awk '{ print  $0, ($6/$3), ($7/$4), ($8/$5)}' | sed 's/ /\t/g' | awk -v alScore=$alScore  '($9>=alScore && $11>=alScore)' | awk '($9>=alScore && $11>=alScore)'  >  $workingDir/tocheckrepeats4.temp.txt
cut -f1 $workingDir/tocheckrepeats4.temp.txt > $workingDir/tocheckrepeats5.temp.txt
grep -Fwf $workingDir/tocheckrepeats5.temp.txt  $workingDir/filtered.bst > $workingDir/tocheckrepeats6.temp.txt
sort -k1,1 -k12,12nr  $workingDir/tocheckrepeats6.temp.txt | sort -k1,1 -u  > $workingDir/tocheckrepeats7.temp.txt
awk -F"\t" -v vectorString=$vectorString -f $scriptDir/extractIS.awk $workingDir/tocheckrepeats7.temp.txt  | grep $vectorString > $workingDir/repeats.resultsNoDup.temp.csv;
sort -k1,1 -u $workingDir/repeats.resultsNoDup.temp.csv > $workingDir/repeats.resultsNoDupSingle.temp.csv

cat $workingDir/repeats.resultsNoDup.csv $workingDir/repeats.resultsNoDup.temp.csv | sort -k1,1 -u >   $workingDir/repeats.resultsNoDup.combined.csv;
mv $workingDir/repeats.resultsNoDup.combined.csv $workingDir/repeats.resultsNoDup.csv 

cat $workingDir/repeats.resultsNoDupSingle.csv $workingDir/repeats.resultsNoDupSingle.temp.csv | sort -k1,1 -u > $workingDir/repeats.resultsNoDupSingle.combined.csv
mv  $workingDir/repeats.resultsNoDupSingle.combined.csv  $workingDir/repeats.resultsNoDupSingle.csv
######################################################################################
######################################################################################
  #get only vec-gen IS reads and remove vec-vec IS reads if any
        cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/repeats.resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2!=vectorString || $3!=vectorString )' > $workingDir/repeats.vecGen-resultsNoDupSingle.csv
        cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/repeats.vecGen-resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2!=vectorString )' > $workingDir/repeats.vecGen-resultsNoDupSingle0.csv

 	awk '{print ($1,$2,$3,$4,$5,$7,$6,$8,$9,$10,$11) }'  $workingDir/repeats.vecGen-resultsNoDupSingle0.csv   >  $workingDir/repeats.vecGen-resultsNoDupSingle0.csv5
        cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/repeats.vecGen-resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2==vectorString )' > $workingDir/repeats.vecGen-resultsNoDupSingle00.csv
         awk '{print ($1,$3,$2,$5,$4,$6,$7,$8,$9,$10,$11) }'   $workingDir/repeats.vecGen-resultsNoDupSingle00.csv >  $workingDir/repeats.vecGen-resultsNoDupSingle9.csv

###############

       echo -e  "SeqID\tChr\tVector\tGenomic_IS\tVector_IS\tGenomic_Strand\tVector_Strand\tSoftClip_Span\tSequence" > $workingDir/repeats.ResultsCompleteUnClustered.csv
        cat  $workingDir/repeats.vecGen-resultsNoDupSingle0.csv5  $workingDir/repeats.vecGen-resultsNoDupSingle9.csv | sed 's/ /\t/g' | cut -f1,2,3,4,5,6,7,9,10,11 >>  $workingDir/repeats.ResultsCompleteUnClustered.csv
sed -i 's/\t/,/g' $workingDir/repeats.ResultsCompleteUnClustered.csv
        cut -d " " -f 1,2,3,4,5,6,7 $workingDir/repeats.resultsNoDupSingle.csv |awk -v vectorName=$vectorString -f $scriptDir/formatIS.awk |sort -k4nr >$workingDir/repeats.resultsNoDup.csv.total

        sort -k1,1 -k2,2n $workingDir/repeats.resultsNoDup.csv.total| awk -f $scriptDir/solveIS.awk > $workingDir/repeats.resultsNoDup.csv.total.results 
#####################################################################################################
#####################################################################################################
#####################################################################################################
cut -d ',' -f1 $workingDir/repeats.ResultsCompleteUnClustered.csv > $workingDir/id.temp
echo -e  "SeqID\tChr\tVector\tGenomic_IS\tVector_IS\tGenomic_Strand\tVector_Strand\tSoftClip_Span\tSequence" > $workingDir/ResultsCompleteUnClustered.csv.1
sed -i 's/\t/,/g'  $workingDir/ResultsCompleteUnClustered.csv.1
grep -v -Fwf $workingDir/id.temp $workingDir/ResultsCompleteUnClustered.csv >> $workingDir/ResultsCompleteUnClustered.csv.1
mv $workingDir/ResultsCompleteUnClustered.csv.1 $workingDir/ResultsCompleteUnClustered.csv

cut -d ',' -f1 $workingDir/repeats.resultsNoDupSingle.csv > $workingDir/id.temp1
grep -v -Fwf $workingDir/id.temp1 $workingDir/resultsNoDupSingle.csv > $workingDir/resultsNoDupSingle.csv.1
mv $workingDir/resultsNoDupSingle.csv.1 $workingDir/resultsNoDupSingle.csv
cut -d " " -f 1,2,3,4,5,6,7 $workingDir/resultsNoDupSingle.csv |awk -v vectorName=$vectorString -f $scriptDir/formatIS.awk |sort -k4nr >$workingDir/resultsNoDup.csv.total
#sort -k1,1 -k2,2n $workingDir/resultsNoDup.csv.total| awk  -f $scriptDir/solveIS.awk > $workingDir/resultsNoDup.csv.total.results
sort -k1,1 -k2,2n $workingDir/resultsNoDup.csv.total| awk -v range=$range -f $scriptDir/solveIS.awk > $workingDir/resultsNoDup.csv.total.results

rm $workingDir/tocheckrepeats* $workingDir/repeats.vecGen-resultsNoDupSingle* $workingDir/filtered.temp.fa $workingDir/filtered.sort.temp.bst
rm $workingDir/repeats.resultsNoDupSingle.temp.csv $workingDir/id.temp $workingDir/id.temp1 $workingDir/repeats.resultsNoDup.temp.csv $workingDir/selected $workingDir/selected1
#####################################################################################################
