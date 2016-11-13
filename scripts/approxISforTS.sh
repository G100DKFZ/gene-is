outputDir=$1
scriptDir=$2
alignmentSam=$3
alignmentNoDup=$4
vectorString=$5
refGenomeOnlyIndexBWA=$6
vecOnlyIndexBWA=$7
alignerBWA=$8
samtools=$9
#APPROX IS for paired end TS data

#exract those reads from nondup.sam where the complete read aligns either with genome/vector
#samtools view $alignmentNoDup | wc -l > $outputDir/lines.txt
#samtools view $alignmentNoDup > $outputDir/completAlignment.sorted.sam
$samtools view  $outputDir/completAlignment.sorted.bam | wc -l > $outputDir/lines.txt
$samtools view  $outputDir/completAlignment.sorted.bam >  $outputDir/completAlignment.sorted.sam



python $scriptDir/approxIS_MeMextract.py  $outputDir/lines.txt  $outputDir/completAlignment.sorted.sam >  $outputDir/approx_MM.sam
sed '/^$/d'  $outputDir/approx_MM.sam >  $outputDir/approx_MM1.sam
#exract vec-gen paired reads
awk -v vectorString=$vectorString '($3==vectorString && $7!="=")' $outputDir/approx_MM1.sam > $outputDir/approx_MM2.sam
awk -v vectorString=$vectorString '($3!=vectorString && $7!="=")' $outputDir/approx_MM1.sam > $outputDir/approx_MM3.sam
cut -f 1 $outputDir/approx_MM2.sam > $outputDir/approx_MM4.sam

awk 'NR==FNR{tgts[$1]; next} $1 in tgts'  $outputDir/approx_MM4.sam $outputDir/approx_MM3.sam > $outputDir/approx_MM5.sam
cut -f 1 $outputDir/approx_MM5.sam > $outputDir/approx_MM6.sam

awk 'NR==FNR{tgts[$1]; next} $1 in tgts' $outputDir/approx_MM6.sam $outputDir/approx_MM2.sam > $outputDir/approx_MM7.sam

cat $outputDir/approx_MM5.sam $outputDir/approx_MM7.sam | sort > $outputDir/approx_MM8.sam
#create fastq
cat  $outputDir/approx_MM7.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' >  $outputDir/approx_MM9.fastq
cat  $outputDir/approx_MM5.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' >  $outputDir/approx_MM10.fastq


$alignerBWA mem -M  $refGenomeOnlyIndexBWA   $outputDir/approx_MM9.fastq  >  $outputDir/approx_MM11.sam
sed -i "/\@SQ\b/d"    $outputDir/approx_MM11.sam

awk '($3=="*")'  $outputDir/approx_MM11.sam | cut -f1  >  $outputDir/approx_MM12.sam
awk '$6  ~ /S/'  $outputDir/approx_MM11.sam | cut -f1 >>  $outputDir/approx_MM12.sam
sort -k1,1 -u  $outputDir/approx_MM12.sam >  $outputDir/approx_MM13.sam

#align gen reads with vec only
$alignerBWA mem -M $vecOnlyIndexBWA  $outputDir/approx_MM10.fastq   > $outputDir/approx_MM14.sam
sed -i "/\@SQ\b/d"   $outputDir/approx_MM14.sam 
awk '($3=="*")' $outputDir/approx_MM14.sam | cut -f1 > $outputDir/approx_MM15.sam
awk '$6  ~ /S/' $outputDir/approx_MM14.sam | cut -f1 >> $outputDir/approx_MM15.sam
sort -k1,1 -u $outputDir/approx_MM15.sam > $outputDir/approx_MM16.sam

awk 'NR==FNR{tgts[$1]; next} $1 in tgts' $outputDir/approx_MM16.sam $outputDir/approx_MM13.sam > $outputDir/approx_MM17.sam
awk 'NR==FNR{tgts[$1]; next} $1 in tgts' $outputDir/approx_MM17.sam $outputDir/approx_MM8.sam > $outputDir/approx_MM18.sam

#########################
cut -f 1 $outputDir/approx_MM18.sam > $outputDir/approx_MM18a.sam
awk 'NR==FNR{tgts[$1]; next} $1 in tgts' $outputDir/approx_MM18a.sam $alignmentSam > $outputDir/approx_MM18b.sam
cut -f1,2,3,4,5,6,7,8,9,10,11,12 $outputDir/approx_MM18b.sam > $outputDir/approx_MM18c.sam
awk -F'\t' -vOFS='\t' '{ gsub("M", "", $6) ; print }'  $outputDir/approx_MM18c.sam > $outputDir/approx_MM18d.sam
less $outputDir/approx_MM18b.sam | wc -l > $outputDir/lines.txt
#get lines to modify by pyhton for addition 
sed -n '1~2p' $outputDir/approx_MM18d.sam > $outputDir/approx_MM18e.sam
#no addition required in these lines
sed -n '0~2p' $outputDir/approx_MM18d.sam > $outputDir/approx_MM18f.sam 

python $scriptDir/approxIS_IScorrection.py  $outputDir/lines.txt  $outputDir/approx_MM18e.sam >  $outputDir/approx_MM18g.sam
paste  $outputDir/approx_MM18e.sam  $outputDir/approx_MM18g.sam >  $outputDir/approx_MM18h.sam

#just to copy column 3 of pos to a col 13 to make it compatibe with corrected IS file
cut -f4  $outputDir/approx_MM18f.sam  >   $outputDir/approx_MM18i.sam 
paste  $outputDir/approx_MM18f.sam   $outputDir/approx_MM18i.sam >   $outputDir/approx_MM18j.sam 


#cat two above files
cat  $outputDir/approx_MM18h.sam  $outputDir/approx_MM18j.sam | sort -k1,1 >  $outputDir/approx_MM18k.sam 
#########################
cut -f 1,2,3,4,6,7,8,10,12,13 $outputDir/approx_MM18k.sam > $outputDir/approx_MM19.sam
sed  -i 's/NM:i://g' $outputDir/approx_MM19.sam
awk -F'\t' -vOFS='\t' '{ gsub("M", "", $5) ; print }' $outputDir/approx_MM19.sam > $outputDir/approx_MM20.sam
awk '{print ((($5- $9)/$5))*100}' $outputDir/approx_MM20.sam > $outputDir/approx_MM21.sam
paste $outputDir/approx_MM20.sam $outputDir/approx_MM21.sam > $outputDir/approx_MM22.sam


awk  -v vectorString=$vectorString  '($3!=vectorString)' $outputDir/approx_MM22.sam | cut -f1,2,3,4,5,8,10,11 | sort -k1,1 > $outputDir/approx_MM23.sam
awk  -v vectorString=$vectorString  '($3==vectorString)' $outputDir/approx_MM22.sam | cut -f1,2,3,4,5,8,10,11 | sort -k1,1 > $outputDir/approx_MM24.sam
paste $outputDir/approx_MM23.sam $outputDir/approx_MM24.sam > $outputDir/approx_MM25.sam


cut -f 1,2,3,4,5,6,7,8,10,11,12,13,14,15,16 $outputDir/approx_MM25.sam > $outputDir/approx_MM26.sam
awk '($8>=95 && $15>=95)'  $outputDir/approx_MM26.sam > $outputDir/approx_MM27.sam

#replace sam flag with strand info.
#sam flag gen 2nd column, vec sam flag is 9th column
 awk -F'\t' -vOFS='\t' '{ gsub("113", "-", $9) ; gsub("177", "-", $9) ; gsub("145", "-", $9) ;gsub("81", "-", $9) ; gsub("65", "+", $9) ; gsub("129", "+", $9) ; gsub("161", "+", $9) ; gsub("97", "+", $9) ; print }' $outputDir/approx_MM27.sam > $outputDir/approx_MM28.sam

 awk -F'\t' -vOFS='\t' '{ gsub("113", "-", $2) ; gsub("177", "-", $2) ; gsub("145", "-", $2) ;gsub("81", "-", $2) ; gsub("65", "+", $2) ; gsub("129", "+", $2) ; gsub("161", "+", $2) ; gsub("97", "+", $2) ; print }' $outputDir/approx_MM28.sam > $outputDir/approx_MM29.sam

#sorted according to sequneces both 
 sort -k6,6 -k13,13 -u $outputDir/approx_MM29.sam > $outputDir/approx_MM30.sam
#sorted according to IS both gen and vector
#Corrected IS are 7 and 14 column so do not use 4 and 11.
sort -k3,3 -k7,7 -k14,14 -u $outputDir/approx_MM30.sam > $outputDir/approx_MM31.sam

echo -e  "SeqID\tChr\tVector\tGenomic_IS\tVector_IS\tGenomic_Strand\tVector_Strand\tGenomeSpan\tVectorSpan\tGenomePercentageIdentity\tVectorPercentageIdentity\tGenomeSeq\tVectorSequence" > $outputDir/Approx.ResultsCompleteUnClustered.csv
awk '{print ($1,$3,$10,$7,$14,$2,$9,$5,$12,$8,$15,$6,$13)}' $outputDir/approx_MM31.sam >> $outputDir/Approx.ResultsCompleteUnClustered.csv
sed -i 's/ /,/g' $outputDir/Approx.ResultsCompleteUnClustered.csv
sed -i 's/\t/,/g' $outputDir/Approx.ResultsCompleteUnClustered.csv
cut -d "," -f1,2,3,4,5,6,7 $outputDir/Approx.ResultsCompleteUnClustered.csv | sed 's/,/ /g' > $outputDir/approx.single.csv
cut -d " " -f 1,2,3,4,5,6,7  $outputDir/approx.single.csv | tail -n+2 |awk -v vectorString=$vectorString -f $scriptDir/approxIS_formatIS.awk |sort -k4nr > $outputDir/approx.single.csv.total
sort -k1,1 -k2,2n $outputDir/approx.single.csv.total | awk -f $scriptDir/solveIS.awk >  $outputDir/approx.resultsNoDup.csv.total.results

rm $outputDir/approx_* $outputDir/lines.txt
#for single end csv require ID, chr2, pos 4th col, 7th col strand
#if chr3, pos5, strand6
##########################
