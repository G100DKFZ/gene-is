alignmentOut=$2
workingDir=$1
scriptDir=$3
refGenome=$4
vectorString=$5
aligner=$6
libDir=$7
refGenomeIndex=$8
minIden=$9
range=${10}
samtools=${11}
specificity=${12}
#### SCLIP extraction

mv $workingDir/completAlignment.nodup.bam   $workingDir/completAlignment.sorted.bam
mv $workingDir/completAlignment.nodup.bam.bai $workingDir/completAlignment.sorted.bam.bai
#echo "perl -I $libDir $scriptDir/extractSClip.pl -i $workingDir/$alignmentOut -o $workingDir/ --min_pct_id 95 --min_pct_hq 90 --ref_genome ${refGenome} ";
#perl -I $libDir $scriptDir/extractSClip.pl -i $workingDir/$alignmentOut -o $workingDir/ --min_pct_id 95 --min_pct_hq 90 --ref_genome ${refGenome};


echo "$samtools view $workingDir/$alignmentOut | awk -v thr=20 -v specific=$specificity -f $scriptDir/extractSClip.awk > $workingDir/$alignmentOut.sclip.txt"
$samtools view $workingDir/$alignmentOut | awk -v thr=20 -v specific=$specificity -f $scriptDir/extractSClip.awk > $workingDir/$alignmentOut.sclip.txt

#### if the aliger selected is blast then each single chromosome is analyzed separatelly for speed
if [ "$aligner" = "blast" ]; then
	for i in {1..22}; do
		chromosomeList[$i]=chr$i;
	done
	chromosomeList[23]=chrX
	chromosomeList[24]=chrY
	chromosomeList[25]=chrM
	chromosomeList[26]=$vectorString
	for n in {1..26}; do
		echo "Processing ${chromosomeList[${n}]}"

                echo "awk -v fa=true -v nodup=true -v chr=${chromosomeList[${n}]} -v outDir=$workingDir -f $scriptDir/filterSCliped.awk $workingDir/$alignmentOut.sclip.txt;"
                awk -v fa=true -v nodup=true -v chr=${chromosomeList[${n}]} -v outDir=$workingDir -f $scriptDir/filterSCliped.awk $workingDir/$alignmentOut.sclip.txt;
                
		echo "megablast -d $refGenomeIndex -i $workingDir/${chromosomeList[${n}]}.fa -m 8 -e 1e-5 > $workingDir/${chromosomeList[${n}]}.bst 2>/dev/null"
		megablast -d $refGenomeIndex -i $workingDir/${chromosomeList[${n}]}.fa -m 8 -e 1e-5 > $workingDir/${chromosomeList[${n}]}.bst 2>/dev/null
                
		echo "awk -F\"\t\" -v vectorStr=$vectorString -f $scriptDir/extractIS.awk $workingDir/${chromosomeList[${n}]}.bst |grep $vectorString >> $workingDir/resultsNoDup.csv"
		awk -F"\t" -v vectorStr=$vectorString -f $scriptDir/extractIS.awk $workingDir/${chromosomeList[${n}]}.bst |grep $vectorString >> $workingDir/resultsNoDup.csv;
	done
# by default blat is used with -m 8 blast tabular option
else
	echo "Processing whole alignment"
	echo "sort -k1,1 -k2,2n -k5,5 $workingDir/$alignmentOut.sclip.txt > $workingDir/$alignmentOut.sclip.txt.app; mv $workingDir/$alignmentOut.sclip.txt.app $workingDir/$alignmentOut.sclip.txt"
	sort -k1,1 -k2,2n -k5,5 $workingDir/$alignmentOut.sclip.txt > $workingDir/$alignmentOut.sclip.txt.app
	mv $workingDir/$alignmentOut.sclip.txt.app $workingDir/$alignmentOut.sclip.txt

##        echo "awk -v fa=true -v nodup=true -v outDir=$workingDir -f $scriptDir/filterSCliped.awk $workingDir/$alignmentOut.sclip.txt;"
# #       awk -v fa=true -v nodup=true -v outDir=$workingDir -f $scriptDir/filterSCliped.awk $workingDir/$alignmentOut.sclip.txt;
 
	sort -k7,7 -u $workingDir/$alignmentOut.sclip.txt > $workingDir/$alignmentOut.sclip.txt1
	mv $workingDir/$alignmentOut.sclip.txt1 $workingDir/$alignmentOut.sclip.txt
	awk '{print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"@"$9"\n"$5 }' $workingDir/$alignmentOut.sclip.txt  > $workingDir/filtered.fa
	        
       
        echo " "
	echo "BLAT alignment in process... "        
        echo "$aligner $refGenomeIndex $workingDir/filtered.fa -out=blast8 -minIdentity=$minIden $workingDir/filtered.bst"
        $aligner $refGenomeIndex $workingDir/filtered.fa -out=blast8 -minIdentity=$minIden $workingDir/filtered.bst

	echo "sort -k1,1 -k12,12nr $workingDir/filtered.bst > $workingDir/filtered.sort.bst"
	sort -k1,1 -k12,12nr $workingDir/filtered.bst > $workingDir/filtered.sort.bst

	echo "awk -F\"\t\" -v vectorStr=$vectorString -f $scriptDir/extractIS.awk $workingDir/filtered.sort.bst |grep $vectorString > $workingDir/resultsNoDup.csv"
	awk -F"\t" -v vectorStr=$vectorString -f $scriptDir/extractIS.awk $workingDir/filtered.sort.bst |grep $vectorString > $workingDir/resultsNoDup.csv;

	echo "sort -k1,1 -u $workingDir/resultsNoDup.csv > $workingDir/resultsNoDupSingle.csv"
	sort -k1,1 -u $workingDir/resultsNoDup.csv > $workingDir/resultsNoDupSingle.csv

	#get only vec-gen IS reads and remove vec-vec IS reads if any
	cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2!=vectorString || $3!=vectorString )' > $workingDir/vecGen-resultsNoDupSingle.csv
	cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/vecGen-resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2!=vectorString )' > $workingDir/vecGen-resultsNoDupSingle0.csv
	

	awk '{print ($1,$2,$3,$4,$5,$7,$6,$8,$9,$10,$11) }'  $workingDir/vecGen-resultsNoDupSingle0.csv   >  $workingDir/vecGen-resultsNoDupSingle0.csv5

        cut -d " "   -f 1,2,3,4,5,6,7,8,9,10-  $workingDir/vecGen-resultsNoDupSingle.csv | awk -v vectorString=$vectorString '($2==vectorString )' > $workingDir/vecGen-resultsNoDupSingle00.csv
         awk '{print ($1,$3,$2,$5,$4,$6,$7,$8,$9,$10,$11) }'   $workingDir/vecGen-resultsNoDupSingle00.csv >  $workingDir/vecGen-resultsNoDupSingle9.csv


	echo -e  "SeqID\tChr\tVector\tGenomic_IS\tVector_IS\tGenomic_Strand\tVector_Strand\tSoftClip_Span\tSequence" > $workingDir/ResultsCompleteUnClustered.csv
	cat  $workingDir/vecGen-resultsNoDupSingle0.csv5  $workingDir/vecGen-resultsNoDupSingle9.csv | sed 's/ /\t/g' | cut -f1,2,3,4,5,6,7,9,10,11 >>  $workingDir/ResultsCompleteUnClustered.csv

sed -i 's/\t/,/g' $workingDir/ResultsCompleteUnClustered.csv

	rm $workingDir/vecGen-resultsNoDupSingle*
	echo " "
	echo "Clustering in process..."
	echo "cut -d \" \" -f 1,2,3,4,5,6,7 $workingDir/resultsNoDupSingle.csv |awk -v vectorName=$vectorString -f $scriptDir/formatIS.awk |sort -k4nr >$workingDir/resultsNoDup.csv.total"
	cut -d " " -f 1,2,3,4,5,6,7 $workingDir/resultsNoDupSingle.csv |awk -v vectorName=$vectorString -f $scriptDir/formatIS.awk |sort -k4nr >$workingDir/resultsNoDup.csv.total

	echo "sort -k1,1 -k2,2n $workingDir/resultsNoDup.csv.total| awk -v range=$range  -f $scriptDir/solveIS.awk > $workingDir/resultsNoDup.csv.total.results"
	sort -k1,1 -k2,2n $workingDir/resultsNoDup.csv.total| awk -v range=$range -f $scriptDir/solveIS.awk > $workingDir/resultsNoDup.csv.total.results 
fi
