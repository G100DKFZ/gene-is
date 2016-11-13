vectorString=$1
blatIndexFiles=$2
workingDir=$3
scriptDir=$4
blatAligner=$5

awk -v vectorString=$vectorString '($2!=vectorString  || $3!=vectorString )' $workingDir/resultsNoDupSingle.csv | cut -d " " -f10,11 | sed 's/ /@/g' > $workingDir/tocheckfalsepositives.txt1
awk -v vectorString=$vectorString '($2!=vectorString  || $3!=vectorString )' $workingDir/resultsNoDupSingle.csv | cut -d " " -f1,3,5 | sed 's/ /@/g' > $workingDir/tocheckfalsepositives.txt2

grep -Fwf $workingDir/tocheckfalsepositives.txt1 $workingDir/filtered.fa >  $workingDir/tocheckfalsepositives.txt3
grep -Fwf $workingDir/tocheckfalsepositives.txt2 $workingDir/tocheckfalsepositives.txt3  > $workingDir/tocheckfalsepositives.txt4

grep -A1 -Fwf  $workingDir/tocheckfalsepositives.txt4 $workingDir/filtered.fa | sed 's/--//g' | sed '/^$/d' >  $workingDir/tocheckfalsepositives.txt5

grep '>' $workingDir/tocheckfalsepositives.txt5 > $workingDir/tocheckfalsepositives.txt6
grep -v '>' $workingDir/tocheckfalsepositives.txt5 > $workingDir/tocheckfalsepositives.txt7

paste  $workingDir/tocheckfalsepositives.txt6   $workingDir/tocheckfalsepositives.txt7 >  $workingDir/tocheckfalsepositives.txt8
cut -d "@" -f 6  $workingDir/tocheckfalsepositives.txt8 >  $workingDir/tocheckfalsepositives.txt9
paste $workingDir/tocheckfalsepositives.txt6 $workingDir/tocheckfalsepositives.txt9 $workingDir/tocheckfalsepositives.txt7 > $workingDir/tocheckfalsepositivestxt10
cat  $workingDir/tocheckfalsepositivestxt10 | wc -l  > $workingDir/Lines 

python $scriptDir/removeSpartFromSequnce.py  $workingDir/Lines $workingDir/tocheckfalsepositivestxt10 > $workingDir/tocheckfalsepositives.txt11.fa


#######
#separate vec and chr ids header sequences as they will be processed separately
grep -A1 $vectorString $workingDir/tocheckfalsepositives.txt11.fa | sed 's/--//g' | sed '/^$/d' > $workingDir/tocheckfalsepositives.txt11a.fa
$blatAligner $blatIndexFiles $workingDir/tocheckfalsepositives.txt11a.fa -out=blast8 -minIdentity=95 $workingDir/tocheckfalsepositives.txt11a.bst
sort -k1,1 -u   $workingDir/tocheckfalsepositives.txt11a.bst | grep -v 'chr' | cut -d "@" -f1 | sort | uniq >  $workingDir/tocheckfalsepositives.txt11a.bst.ed1
# seprate chr IDS
grep -A1 '@chr'  $workingDir/tocheckfalsepositives.txt11.fa | sed 's/--//g' | sed '/^$/d' > $workingDir/tocheckfalsepositives.txt11b.fa
$blatAligner $blatIndexFiles $workingDir/tocheckfalsepositives.txt11b.fa -out=blast8 -minIdentity=95 $workingDir/tocheckfalsepositives.txt11b.bst
#then grep above ids from resultsNoDupSingle.csv but this also results in adding some seq which should not be there after extra file.
#So now I comment above and dont do this change but cna keep in mind that this can sometime result in losing some IS also
#so I can rescue some IS still if these ids were in filtered.bst only once then rscue them..may be improve later
cut -f1,3,4  $workingDir/tocheckfalsepositives.txt11b.bst | uniq -c  | sed 's/^ *//g' | sed 's/ /\t/' >  $workingDir/tocheckfalsepositives.txt11b.bst.ed1
cut -f1  $workingDir/tocheckfalsepositives.txt11b.bst.ed1 >  $workingDir/tocheckfalsepositives.txt11b.bst.ed2
cut -f2,3,4  $workingDir/tocheckfalsepositives.txt11b.bst.ed1 >  $workingDir/tocheckfalsepositives.txt11b.bst.ed3
 paste  $workingDir/tocheckfalsepositives.txt11b.bst.ed3  $workingDir/tocheckfalsepositives.txt11b.bst.ed2 >  $workingDir/tocheckfalsepositives.txt11b.bst.ed4
sort -k1,1 -u  $workingDir/tocheckfalsepositives.txt11b.bst.ed4 |  awk '($4==1)' | cut -d "@" -f1 | sort | uniq  >  $workingDir/tocheckfalsepositives.txt11b.bst.ed5

#combined correct all
cut -d "@"  -f1  $workingDir/tocheckfalsepositives.txt2 >  $workingDir/tocheckfalsepositives.txt1.ids.ed
cat  $workingDir/tocheckfalsepositives.txt11a.bst  $workingDir/tocheckfalsepositives.txt11b.bst | cut -d "@" -f1 | sort | uniq >  $workingDir/tocheckfalsepositives.txt1ab.bst.ids
grep -v -Fwf  $workingDir/tocheckfalsepositives.txt1ab.bst.ids  $workingDir/tocheckfalsepositives.txt1.ids.ed >  $workingDir/tocheckfalsepositives.txt.rescueIds

cat  $workingDir/tocheckfalsepositives.txt11b.bst.ed5  $workingDir/tocheckfalsepositives.txt11a.bst.ed1 $workingDir/tocheckfalsepositives.txt.rescueIds  >  $workingDir/tocheckfalsepositives.txt11ab.ed

# generate results files again after extra filtering
grep -Fwf $workingDir/tocheckfalsepositives.txt11ab.ed $workingDir/resultsNoDupSingle.csv  > $workingDir/resultsNoDupSingle.csv.extraFiltered
mv $workingDir/resultsNoDupSingle.csv  $workingDir/resultsNoDupSingle.csv.full
mv $workingDir/resultsNoDupSingle.csv.extraFiltered $workingDir/resultsNoDupSingle.csv  

echo -e  "SeqID\tChr\tVector\tGenomic_IS\tVector_IS\tGenomic_Strand\tVector_Strand\tSoftClip_Span\tSequence" > $workingDir/ResultsCompleteUnClustered.extraFiltered.csv
sed -i 's/,/\t/g' $workingDir/ResultsCompleteUnClustered.csv 
grep -Fwf  $workingDir/tocheckfalsepositives.txt11ab.ed   $workingDir/ResultsCompleteUnClustered.csv >>  $workingDir/ResultsCompleteUnClustered.extraFiltered.csv
mv $workingDir/ResultsCompleteUnClustered.extraFiltered.csv $workingDir/ResultsCompleteUnClustered.csv
sed -i 's/\t/,/g'  $workingDir/ResultsCompleteUnClustered.csv
cut -d " " -f 1,2,3,4,5,6,7 $workingDir/resultsNoDupSingle.csv |awk -v vectorName=$vectorString -f $scriptDir/formatIS.awk |sort -k4nr >$workingDir/resultsNoDup.csv.total
sort -k1,1 -k2,2n $workingDir/resultsNoDup.csv.total| awk -f $scriptDir/solveIS.awk > $workingDir/resultsNoDup.csv.total.results
rm $workingDir/resultsNoDupSingle.csv.full $workingDir/tocheckfalsepositives* $workingDir/Lines*
#################
