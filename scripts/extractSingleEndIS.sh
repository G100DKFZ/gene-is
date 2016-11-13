#!/bin/bash
#################################################################
 ## To process single end alignment file for processing IS
################################################################
alignmentOut=$1 
workingDir=$2 
scriptDir=$3 
vectorString=$4
vectorIndexFiles=$5
referenceIndexFiles=$6
forwardOut=$7
blatAligner=$8
aligner=$9
refGenomeIndex=${10}
minIden=${11}
range=${12}
#################################################################
#process file to get reads that align with vector and genome
echo "sed -i "/\@SQ\b/d"  $workingDir/$alignmentOut"
sed -i "/\@SQ\b/d"  $workingDir/$alignmentOut
#from complete alignment extract those SE reads that have ¿S¿part

echo "awk '$6 ~ /S/'  $workingDir/$alignmentOut >   $workingDir/completAlignment_onlyS.sam"
awk '$6 ~ /S/'  $workingDir/$alignmentOut >   $workingDir/completAlignment_onlyS.sam

echo "sed -i '/^$/d'  $workingDir/completAlignment_onlyS.sam"
sed -i '/^$/d'  $workingDir/completAlignment_onlyS.sam



#remove ambiguous CIGAR raeds
echo "CountLines=$(sed -n '$=' $workingDir/completAlignment_onlyS.sam)"
CountLines=$(sed -n '$=' $workingDir/completAlignment_onlyS.sam)
echo "echo $CountLines > $workingDir/Lines.txt"
echo $CountLines > $workingDir/Lines.txt

echo "python $scriptDir/CIGAR_MS_Correction_Oct2014.py $workingDir/Lines.txt $workingDir/completAlignment_onlyS.sam > $workingDir/completAlignment_onlyS_corrected.sam"
python $scriptDir/CIGAR_MS_Correction_Oct2014.py $workingDir/Lines.txt $workingDir/completAlignment_onlyS.sam > $workingDir/completAlignment_onlyS_corrected.sam 
#get unique reads

echo "sort -u -k1,1  $workingDir/completAlignment_onlyS_corrected.sam >  $workingDir/completAlignment_onlyS_corrected_sorted.sam"
sort -u -k1,1  $workingDir/completAlignment_onlyS_corrected.sam >  $workingDir/completAlignment_onlyS_corrected_sorted.sam
#extract desired columns

echo "cut -f 1,2,3,4,5,6,7,8,9,10,11,12 $workingDir/completAlignment_onlyS_corrected_sorted.sam > $workingDir/completAlignment_onlyS_corrected_sorted1.sam"
cut -f 1,2,3,4,5,6,7,8,9,10,11,12 $workingDir/completAlignment_onlyS_corrected_sorted.sam > $workingDir/completAlignment_onlyS_corrected_sorted1.sam
# to get strand information replace sam flag with related strand

echo "awk '{ sub(/0$/, "+", $2) }1' $workingDir/completAlignment_onlyS_corrected_sorted1.sam | awk '{ sub(/256$/, "+", $2) }1'  | awk '{ sub(/16$/, "-", $2) }1' | awk '{ sub(/272$/, "-", $2) }1'  |  awk '{ sub(/0$/, "+", $7) }1'   | awk '{ sub(/256$/, "+", $7) }1'  | awk '{ sub(/16$/, "-", $7) }1' | awk '{ sub(/272$/, "-", $7) }1'  > $workingDir/completAlignment_onlyS_corrected_sorted1_strand.sam"

awk '{ sub(/0$/, "+", $2) }1' $workingDir/completAlignment_onlyS_corrected_sorted1.sam | awk '{ sub(/256$/, "+", $2) }1'  | awk '{ sub(/16$/, "-", $2) }1' | awk '{ sub(/272$/, "-", $2) }1'  |  awk '{ sub(/0$/, "+", $7) }1'   | awk '{ sub(/256$/, "+", $7) }1'  | awk '{ sub(/16$/, "-", $7) }1' | awk '{ sub(/272$/, "-", $7) }1'  > $workingDir/completAlignment_onlyS_corrected_sorted1_strand.sam

echo "awk -v OFS="\t" '$1=$1'  $workingDir/completAlignment_onlyS_corrected_sorted1_strand.sam > $workingDir/completAlignment_onlyS_corrected_sorted1_strand1.sam"
awk -v OFS="\t" '$1=$1'  $workingDir/completAlignment_onlyS_corrected_sorted1_strand.sam > $workingDir/completAlignment_onlyS_corrected_sorted1_strand1.sam

#Separate those reads which have CIGAR MS/SM for further processing
echo "python $scriptDir/CIGAR_SeprateSM_Oct2014.py $workingDir/Lines.txt $workingDir/completAlignment_onlyS_corrected_sorted1_strand1.sam > $workingDir/completAlignment_onlyS_SM.sam"
python $scriptDir/CIGAR_SeprateSM_Oct2014.py $workingDir/Lines.txt $workingDir/completAlignment_onlyS_corrected_sorted1_strand1.sam > $workingDir/completAlignment_onlyS_SM.sam

echo "sed -i '/^$/d'  $workingDir/completAlignment_onlyS_SM.sam"
sed -i '/^$/d'  $workingDir/completAlignment_onlyS_SM.sam

echo "python $scriptDir/CIGAR_SeprateMS_Oct2014.py $workingDir/Lines.txt $workingDir/completAlignment_onlyS_corrected_sorted1_strand1.sam >  $workingDir/completAlignment_onlyS_MS.sam"
python $scriptDir/CIGAR_SeprateMS_Oct2014.py $workingDir/Lines.txt $workingDir/completAlignment_onlyS_corrected_sorted1_strand1.sam >  $workingDir/completAlignment_onlyS_MS.sam

echo "sed -i '/^$/d'  $workingDir/completAlignment_onlyS_MS.sam"
sed -i '/^$/d'  $workingDir/completAlignment_onlyS_MS.sam


#extract S part, span and sequence.
echo " python $scriptDir/CIGAR_SM_Sextraction.py   $workingDir/Lines.txt $workingDir/completAlignment_onlyS_SM.sam > $workingDir/completAlignment_onlyS_SM_Sextracted.sam
"
python $scriptDir/CIGAR_SM_Sextraction.py   $workingDir/Lines.txt $workingDir/completAlignment_onlyS_SM.sam > $workingDir/completAlignment_onlyS_SM_Sextracted.sam
echo "python $scriptDir/CIGAR_MS_Sextraction.py  $workingDir/Lines.txt $workingDir/completAlignment_onlyS_MS.sam >  $workingDir/completAlignment_onlyS_MS_Sextracted.sam"
python $scriptDir/CIGAR_MS_Sextraction.py  $workingDir/Lines.txt $workingDir/completAlignment_onlyS_MS.sam >  $workingDir/completAlignment_onlyS_MS_Sextracted.sam
echo " awk -v OFS="\t" '$1=$1'  $workingDir/completAlignment_onlyS_MS_Sextracted.sam >  $workingDir/completAlignment_onlyS_MS_Sextracted1.sam"
awk -v OFS="\t" '$1=$1'  $workingDir/completAlignment_onlyS_MS_Sextracted.sam >  $workingDir/completAlignment_onlyS_MS_Sextracted1.sam

echo "awk -v OFS="\t" '$1=$1' $workingDir/completAlignment_onlyS_SM_Sextracted.sam > $workingDir/completAlignment_onlyS_SM_Sextracted1.sam"
awk -v OFS="\t" '$1=$1' $workingDir/completAlignment_onlyS_SM_Sextracted.sam > $workingDir/completAlignment_onlyS_SM_Sextracted1.sam

#get those reads having length greater than 20
echo "awk '($13 >= 20 )' $workingDir/completAlignment_onlyS_MS_Sextracted1.sam > $workingDir/completAlignment_onlyS_MS_Sextracted1aa.sam"
awk '($13 >= 20 )' $workingDir/completAlignment_onlyS_MS_Sextracted1.sam > $workingDir/completAlignment_onlyS_MS_Sextracted1aa.sam

echo " awk '($13 >= 20 )' $workingDir/completAlignment_onlyS_SM_Sextracted1.sam  > $workingDir/completAlignment_onlyS_SM_Sextracted1aa.sam"
awk '($13 >= 20 )' $workingDir/completAlignment_onlyS_SM_Sextracted1.sam  > $workingDir/completAlignment_onlyS_SM_Sextracted1aa.sam


#####################################################################
#To further refine the results to get rid of all false positives, first align it with VECTOR only sequence. 
#And get all the unique aligned reads and then align with REFERENCE genome only and get uniquely mapped reads.
echo "cut -f1 $workingDir/completAlignment_onlyS_MS_Sextracted1aa.sam >  $workingDir/Exact-ids.txt"
cut -f1 $workingDir/completAlignment_onlyS_MS_Sextracted1aa.sam >  $workingDir/Exact-ids.txt

echo "awk -vExact="$workingDir/Exact-ids.txt" 'BEGIN{while((getline<Exact)>0)l["@"$1]=1}NR%2==1{f=l[$1]?1:0}f' $workingDir/$forwardOut > $workingDir/Exact.fastq
$aligner mem -M $vectorIndexFiles $workingDir/Exact.fastq > $workingDir/ExactVecOnly.sam"
awk -vExact="$workingDir/Exact-ids.txt" 'BEGIN{while((getline<Exact)>0)l[">"$1]=1}NR%2==1{f=l[$1]?1:0}f' $workingDir/$forwardOut > $workingDir/Exact.fastq
$aligner mem -M $vectorIndexFiles $workingDir/Exact.fastq  > $workingDir/ExactVecOnly.sam

echo "sort $workingDir/ExactVecOnly.sam | uniq > $workingDir/ExactVecOnly1.sam"
sort $workingDir/ExactVecOnly.sam | uniq > $workingDir/ExactVecOnly1.sam

echo "cut -f1 $workingDir/ExactVecOnly1.sam | sort | uniq -u > $workingDir/ExactVecOnly.sam.ids"
cut -f1 $workingDir/ExactVecOnly1.sam | sort | uniq -u > $workingDir/ExactVecOnly.sam.ids

echo "awk -vExactVecOnly="$workingDir/ExactVecOnly.sam.ids"  'BEGIN{while((getline<ExactVecOnly)>0)l["@"$1]=1}NR%2==1{f=l[$1]?1:0}f'   $workingDir/$forwardOut > $workingDir/Exact1.fastq"
awk -vExactVecOnly="$workingDir/ExactVecOnly.sam.ids"  'BEGIN{while((getline<ExactVecOnly)>0)l[">"$1]=1}NR%2==1{f=l[$1]?1:0}f'   $workingDir/$forwardOut > $workingDir/Exact1.fastq

echo "$aligner mem -M $referenceIndexFiles $workingDir/Exact1.fastq  > $workingDir/ExactGenOnly.sam"
$aligner mem -M $referenceIndexFiles $workingDir/Exact1.fastq   > $workingDir/ExactGenOnly.sam

echo "sort $workingDir/ExactGenOnly.sam | uniq > $workingDir/ExactGenOnly1.sam"
sort $workingDir/ExactGenOnly.sam | uniq > $workingDir/ExactGenOnly1.sam

echo "cut -f1 $workingDir/ExactGenOnly1.sam | sort | uniq -u  > $workingDir/ExactGenOnly.sam.ids"
cut -f1 $workingDir/ExactGenOnly1.sam | sort | uniq -u  > $workingDir/ExactGenOnly.sam.ids
echo "awk 'NR==FNR{tgts[$1]; next} $1 in tgts' $workingDir/ExactGenOnly.sam.ids $workingDir/completAlignment_onlyS_MS_Sextracted1aa.sam >   $workingDir/completAlignment_onlyS_MS_Sextracted1a.sam"
awk 'NR==FNR{tgts[$1]; next} $1 in tgts' $workingDir/ExactGenOnly.sam.ids $workingDir/completAlignment_onlyS_MS_Sextracted1aa.sam >   $workingDir/completAlignment_onlyS_MS_Sextracted1a.sam

############
echo "cut -f1 $workingDir/completAlignment_onlyS_SM_Sextracted1aa.sam  >  $workingDir/Exact-ids.txt"
cut -f1 $workingDir/completAlignment_onlyS_SM_Sextracted1aa.sam  >  $workingDir/Exact-ids.txt

echo "awk -vExact="$workingDir/Exact-ids.txt" 'BEGIN{while((getline<Exact)>0)l["@"$1]=1}NR%4==1{f=l[$1]?1:0}f' $workingDir/$forwardOut > $workingDir/Exact.fastq
$aligner mem -M $vectorIndexFiles $workingDir/Exact.fastq  > $workingDir/ExactVecOnly.sam"

awk -vExact="$workingDir/Exact-ids.txt" 'BEGIN{while((getline<Exact)>0)l["@"$1]=1}NR%2==1{f=l[$1]?1:0}f' $workingDir/$forwardOut > $workingDir/Exact.fastq
$aligner mem -M $vectorIndexFiles $workingDir/Exact.fastq  > $workingDir/ExactVecOnly.sam

echo "sort $workingDir/ExactVecOnly.sam | uniq > $workingDir/ExactVecOnly1.sam"
sort $workingDir/ExactVecOnly.sam | uniq > $workingDir/ExactVecOnly1.sam
echo "cut -f1 $workingDir/ExactVecOnly1.sam | sort | uniq -u > $workingDir/ExactVecOnly.sam.ids"
cut -f1 $workingDir/ExactVecOnly1.sam | sort | uniq -u > $workingDir/ExactVecOnly.sam.ids

echo "awk -vExactVecOnly="$workingDir/ExactVecOnly.sam.ids"  'BEGIN{while((getline<ExactVecOnly)>0)l["@"$1]=1}NR%4==1{f=l[$1]?1:0}f'   $workingDir/$forwardOut > $workingDir/Exact1.fastq"
awk -vExactVecOnly="$workingDir/ExactVecOnly.sam.ids"  'BEGIN{while((getline<ExactVecOnly)>0)l["@"$1]=1}NR%2==1{f=l[$1]?1:0}f'   $workingDir/$forwardOut > $workingDir/Exact1.fastq

echo "$aligner mem -M $referenceIndexFiles $workingDir/Exact1.fastq > $workingDir/ExactGenOnly.sam"
$aligner mem -M $referenceIndexFiles $workingDir/Exact1.fastq  > $workingDir/ExactGenOnly.sam

echo "sort $workingDir/ExactGenOnly.sam | uniq > $workingDir/ExactGenOnly1.sam"
sort $workingDir/ExactGenOnly.sam | uniq > $workingDir/ExactGenOnly1.sam

echo "cut -f1 $workingDir/ExactGenOnly1.sam | sort | uniq -u  > $workingDir/ExactGenOnly.sam.ids"
cut -f1 $workingDir/ExactGenOnly1.sam | sort | uniq -u  > $workingDir/ExactGenOnly.sam.ids

echo "awk 'NR==FNR{tgts[$1]; next} $1 in tgts' $workingDir/ExactGenOnly.sam.ids $workingDir/completAlignment_onlyS_SM_Sextracted1aa.sam   >  $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam"
awk 'NR==FNR{tgts[$1]; next} $1 in tgts' $workingDir/ExactGenOnly.sam.ids $workingDir/completAlignment_onlyS_SM_Sextracted1aa.sam   >  $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam

#####################################################################

# correct IS in MS
echo "python $scriptDir/CIGAR_MS_MPosCorrectionOct2014.py $workingDir/Lines.txt $workingDir/completAlignment_onlyS_MS_Sextracted1a.sam > $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos0.sam"
python $scriptDir/CIGAR_MS_MPosCorrectionOct2014.py $workingDir/Lines.txt $workingDir/completAlignment_onlyS_MS_Sextracted1a.sam > $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos0.sam

#space to tab before column 15
echo "sed -e "s/\ /\t/g"   $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos0.sam > $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam"
sed -e "s/\ /\t/g"   $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos0.sam > $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam

#add column 1 in MS as column 16
echo "awk -F"\t" 'BEGIN { OFS = "\t" } {$16="1"; print}' $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam  > $workingDir/MS_MS.sam"
awk -F"\t" 'BEGIN { OFS = "\t" } {$16="1"; print}' $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam  > $workingDir/MS_MS.sam

cut -f 1,3 $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam > $workingDir/MS_idChr.txt
cut -f 2 $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam > $workingDir/MS_strand.txt
cut -f 13  $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam > $workingDir/MS_Sspan.txt
cut -f 10 $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam > $workingDir/MS_read.txt

cut -f 15 $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam > $workingDir/MS_IS.txt
# add a new column with 1 at the end to ensure it is MS
cut -f 16 $workingDir/MS_MS.sam >  $workingDir/MS_MS.txt

paste $workingDir/MS_idChr.txt $workingDir/MS_IS.txt $workingDir/MS_strand.txt $workingDir/MS_Sspan.txt $workingDir/MS_read.txt $workingDir/MS_MS.txt | sed  's/\t/@/g' >  $workingDir/MS_header.txt
cut -f 14 $workingDir/completAlignment_onlyS_MS_Sextracted1a_correctedMpos.sam > $workingDir/MS_Sseq.txt

paste  $workingDir/MS_header.txt $workingDir/MS_Sseq.txt | sed -e 's/^/>/' |  sed 's/\t/\n/g' >  $workingDir/completAlignment_S_MS.fa



echo "awk -F"\t" 'BEGIN { OFS = "\t" } {$15="0"; print}' $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam  > $workingDir/SM_SM.sam
cut -f 15   $workingDir/SM_SM.sam > $workingDir/SM_SM.txt"
awk -F"\t" 'BEGIN { OFS = "\t" } {$15="0"; print}' $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam  > $workingDir/SM_SM.sam
cut -f 15   $workingDir/SM_SM.sam > $workingDir/SM_SM.txt

cut -f 1,3,4 $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam > $workingDir/SM_idChrIS.txt
cut -f 2 $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam > $workingDir/SM_strand.txt
cut -f 13  $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam > $workingDir/SM_Sspan.txt
cut -f 10  $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam > $workingDir/SM_read.txt

echo "paste $workingDir/SM_idChrIS.txt $workingDir/SM_strand.txt $workingDir/SM_Sspan.txt $workingDir/SM_read.txt $workingDir/SM_SM.txt |  sed  's/\t/@/g' >  $workingDir/SM_header.txt"
paste $workingDir/SM_idChrIS.txt $workingDir/SM_strand.txt $workingDir/SM_Sspan.txt $workingDir/SM_read.txt $workingDir/SM_SM.txt |  sed  's/\t/@/g' >  $workingDir/SM_header.txt

echo "cut -f 14 $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam > $workingDir/SM_Sseq.txt"
cut -f 14 $workingDir/completAlignment_onlyS_SM_Sextracted1a.sam > $workingDir/SM_Sseq.txt

echo "paste $workingDir/SM_header.txt $workingDir/SM_Sseq.txt | sed -e 's/^/>/' |  sed 's/\t/\n/g' >   $workingDir/completAlignment_S_SM.fa"
paste $workingDir/SM_header.txt $workingDir/SM_Sseq.txt | sed -e 's/^/>/' |  sed 's/\t/\n/g' >   $workingDir/completAlignment_S_SM.fa

#combine MS and SM to make one file
echo "cat $workingDir/completAlignment_S_MS.fa    $workingDir/completAlignment_S_SM.fa >  $workingDir/completAlignment_S.fa"
cat $workingDir/completAlignment_S_MS.fa    $workingDir/completAlignment_S_SM.fa >  $workingDir/completAlignment_S.fa  

echo "$blatAligner $refGenomeIndex  $workingDir/completAlignment_S.fa -out=blast8 -minIdentity=$minIden $workingDir/completAlignment_S.bst"
$blatAligner $refGenomeIndex  $workingDir/completAlignment_S.fa -out=blast8 -minIdentity=$minIden $workingDir/completAlignment_S.bst
echo "echo $vectorString"
echo $vectorString
echo "awk -F"\t" -v vectorStr=$vectorString -f $scriptDir/extractIS.awk  $workingDir/completAlignment_S.bst | grep $vectorString > $workingDir/resultsNoDup.csv"
awk -F"\t" -v vectorStr=$vectorString -f $scriptDir/extractIS.awk  $workingDir/completAlignment_S.bst | grep $vectorString > $workingDir/resultsNoDup.csv


mv $workingDir/completAlignment_S.fa $workingDir/filtered.fa
mv $workingDir/completAlignment_S.bst $workingDir/filtered.bst 
sort -k1,1 -k12,12nr $workingDir/filtered.bst > $workingDir/filtered.sort.bst

rm $workingDir/completAlignment_* $workingDir/SM_* $workingDir/MS_* $workingDir/Exact* $workingDir/Lines*

echo "sort -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 -u $workingDir/resultsNoDup.csv > $workingDir/resultsNoDupSingle0.csv"
sort -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10  -u $workingDir/resultsNoDup.csv > $workingDir/resultsNoDupSingle0.csv
echo "sort -k1,1 -u $workingDir/resultsNoDupSingle0.csv > $workingDir/resultsNoDupSingle.csv"
sort -k1,1 -u $workingDir/resultsNoDupSingle0.csv   > $workingDir/resultsNoDupSingle.csv

rm $workingDir/resultsNoDupSingle0.csv

  #get only vec-gen IS reads and remove vec-vec IS reads if any
 cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2!=vectorString || $3!=vectorString )' > $workingDir/vecGen-resultsNoDupSingle.csv
cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/vecGen-resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2!=vectorString )' > $workingDir/vecGen-resultsNoDupSingle0.csv

awk '{print ($1,$2,$3,$4,$5,$7,$6,$8,$9,$10,$11) }'  $workingDir/vecGen-resultsNoDupSingle0.csv   >  $workingDir/vecGen-resultsNoDupSingle0.csv5
cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/vecGen-resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2==vectorString )' > $workingDir/vecGen-resultsNoDupSingle00.csv
awk '{print ($1,$3,$2,$5,$4,$6,$7,$8,$9,$10,$11) }'   $workingDir/vecGen-resultsNoDupSingle00.csv >  $workingDir/vecGen-resultsNoDupSingle9.csv

echo -e  "SeqID\tChr\tVector\tGenomic_IS\tVector_IS\tGenomic_Strand\tVector_Strand\tSoftClip_Span\tSequence" > $workingDir/ResultsCompleteUnClustered.csv
cat  $workingDir/vecGen-resultsNoDupSingle0.csv5  $workingDir/vecGen-resultsNoDupSingle9.csv | sed 's/ /\t/g' | cut -f1,2,3,4,5,6,7,9,10,11 >>  $workingDir/ResultsCompleteUnClustered.csv
sed -i  's/\t/,/g' $workingDir/ResultsCompleteUnClustered.csv

rm $workingDir/vecGen-resultsNoDupSingle* 

echo "cut -d " " -f 1,2,3,4,5,6,7   $workingDir/resultsNoDupSingle.csv | awk -v vectorName=$vectorString -f $scriptDir/formatIS.awk | sort -k4nr >  $workingDir/resultsNoDup.csv.total"
cut -d " " -f 1,2,3,4,5,6,7   $workingDir/resultsNoDupSingle.csv | awk -v vectorName=$vectorString -f $scriptDir/formatIS.awk | sort -k4nr >  $workingDir/resultsNoDup.csv.total
echo "sort -k1,1 -k2,2n $workingDir/resultsNoDup.csv.total | awk -v range=$range  -f $scriptDir/solveIS.awk >  $workingDir/resultsNoDup.csv.total.results"
sort -k1,1 -k2,2n $workingDir/resultsNoDup.csv.total | awk -v range=$range -f $scriptDir/solveIS.awk >  $workingDir/resultsNoDup.csv.total.results

################################################################################################
