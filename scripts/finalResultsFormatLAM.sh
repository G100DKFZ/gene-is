#!/bin/bash
outputDir=$1
sampleName=$2
len=$3
alScore=$4
idClus=$5
range=$6
notMatchThreshold=$7

less  $outputDir/GeneralStatistics.txt > $outputDir/$sampleName.$len-$alScore-$idClus-$range-$notMatchThreshold.GeneralStatistics.txt
less  $outputDir/ResultsClusteredAnnotated.csv > $outputDir/$sampleName.$len-$alScore-$idClus-$range-$notMatchThreshold.ResultsClusteredAnnotated.csv
mv  $outputDir/wrongBC_.fastq   $outputDir/$sampleName.$len-$alScore-$idClus-$range-$notMatchThreshold.wrongBC_.fastq
less  $outputDir/ResultsTenStrongestClones.csv  > $outputDir/$sampleName.$len-$alScore-$idClus-$range-$notMatchThreshold.ResultsTenStrongestClones.csv
less  $outputDir/repeats.ResultsClusteredAnnotated.csv > $outputDir/$sampleName.$len-$alScore-$idClus-$range-$notMatchThreshold.repeats.ResultsClusteredAnnotated.csv


echo "Chr,Genomic_Position,Strand,Genomic_Span,Read_Sequence,Seq_ID,Sample_Name" >  $outputDir/allSort3
grep '>' $outputDir/*_fa | cut -d ">" -f1 | sed 's/://g' | rev |  cut -d "/" -f1 | rev   >  $outputDir/allSort1
grep '>' $outputDir/*_fa | cut -d ">" -f2- > $outputDir/allSort2
paste  $outputDir/allSort2  $outputDir/allSort1 | sed 's/\t/,/g'  | sed 's/ /,/g' | cut -d"," -f1,2,3,4,5,6,8  >>  $outputDir/allSort3

echo "Chr,Genomic_Position,Strand,Genomic_Span,Read_Sequence,Seq_ID,Sample_Name" >  $outputDir/allSort3r
grep '>' $outputDir/*r1_fa.r | cut -d ">" -f1 | sed 's/://g' | rev |  cut -d "/" -f1 | rev   >  $outputDir/allSort1r
grep '>' $outputDir/*r1_fa.r | cut -d ">" -f2- > $outputDir/allSort2r
paste  $outputDir/allSort2r  $outputDir/allSort1r | sed 's/\t/,/g'  | sed 's/ /,/g'   | cut -d"," -f1,2,3,4,5,6,8  >>  $outputDir/allSort3r

less  $outputDir/allSort3  > $outputDir/$sampleName.$len-$alScore-$idClus-$range-$notMatchThreshold.ResultsCompleteUnclustered.csv
less  $outputDir/allSort3r  > $outputDir/$sampleName.$len-$alScore-$idClus-$range-$notMatchThreshold.repeats.ResultsCompleteUnclustered.csv

rm $outputDir/allSort3 $outputDir/allSort1 $outputDir/allSort2 $outputDir/allSort1r  $outputDir/allSort3r $outputDir/allSort2r 
rm $outputDir/*is..r*

rm   $outputDir/resultsNoDup.csv.total.results  $outputDir/resultsNoDup.csv.total.results.r  $outputDir/testt1 $outputDir/testt2 $outputDir/testt3 $outputDir/tt1 $outputDir/tt2 $outputDir/tt3   $outputDir/GeneralStatistics.txt $outputDir/ResultsClusteredAnnotated.csv     $outputDir/ResultsTenStrongestClones.csv $outputDir/result_AnnotatedIS.csv $outputDir/repeats.ResultsClusteredAnnotated.csv 
