#General Statistics
output_dir=$1
forward=$2
len=$3
alScore=$4
idClus=$5
range=$6
notMatchThreshold=$7

echo "##################################################################" >  $output_dir/GeneralStatistics.txt
echo "GENERAL STATISTICS" >>  $output_dir/GeneralStatistics.txt
echo "##################################################################" >>  $output_dir/GeneralStatistics.txt
echo " " >>  $output_dir/GeneralStatistics.txt
echo " " >>  $output_dir/GeneralStatistics.txt
echo "Number of raw read pairs" >>  $output_dir/GeneralStatistics.txt
F=$(zcat $forward | wc -l)
#F=$(sed -n '$=' $forward)
N=4 
DF=$((F/N))
echo $DF >>  $output_dir/GeneralStatistics.txt
echo " " >>  $output_dir/GeneralStatistics.txt
echo "Number of filtered and trimmed read pairs" >>  $output_dir/GeneralStatistics.txt
TF=$(cat $output_dir/*_.fastq.trimmed2 | wc -l )
#TF=$(sed -n '$=' $output_dir/filtTrim-pair1.fastq)
TN=4
TDF=$((TF/TN))
echo $TDF >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt

echo "Number of unclustered Integration Sites (sequence_count) " >> $output_dir/GeneralStatistics.txt
sed 's/^ *//g'  $output_dir/ResultsClusteredAnnotated.csv | tail -n +2 | sed  's/ /\t/g' | cut -d "," -f 4 | awk '{ sum+=$1} END {print sum}' >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt
echo "Number of clustered Integration Sites" >> $output_dir/GeneralStatistics.txt
cat  $output_dir/ResultsClusteredAnnotated.csv | tail -n +2 | wc -l >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt

echo " " >> $output_dir/GeneralStatistics.txt
echo "Number of unclustered MULTIPLE ALIGNED/REPEATS Integration Sites (sequence_count) " >> $output_dir/GeneralStatistics.txt
sed 's/^ *//g'  $output_dir/repeats.ResultsClusteredAnnotated.csv | tail -n +2 | sed  's/ /\t/g' | cut -d "," -f 4 | awk '{ sum+=$1} END {print sum}' >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt
echo "Number of clustered MULTIPLE ALIGNED/REPEATS Integration Sites" >> $output_dir/GeneralStatistics.txt
cat  $output_dir/repeats.ResultsClusteredAnnotated.csv | tail -n +2 | wc -l >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt



echo "#################################################################" >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt
echo "Wrong barcode combinations detected" >> $output_dir/GeneralStatistics.txt
echo "(Ten highest counts of wrong barcode combinations are shown here)" >> $output_dir/GeneralStatistics.txt
echo "(For complete list of wrong barcodes see file 'wrongBC_.fastq')" >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt
echo "BC1             BC2           Count" >> $output_dir/GeneralStatistics.txt
sort -nr -k3,3 $output_dir/wrongBC_.fastq | head >> $output_dir/GeneralStatistics.txt 
echo " " >> $output_dir/GeneralStatistics.txt
echo "#################################################################" >> $output_dir/GeneralStatistics.txt

echo " " >> $output_dir/GeneralStatistics.txt
echo "Parameters used for analysis" >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt
echo "len: " >> $output_dir/GeneralStatistics.txt
echo $len >> $output_dir/GeneralStatistics.txt

echo " " >> $output_dir/GeneralStatistics.txt
echo "alScore:" >> $output_dir/GeneralStatistics.txt
echo $alScore  >> $output_dir/GeneralStatistics.txt

#echo " " >> $output_dir/GeneralStatistics.txt
#echo "idClus:" >> $output_dir/GeneralStatistics.txt
#echo $idClus >> $output_dir/GeneralStatistics.txt

echo " " >> $output_dir/GeneralStatistics.txt
echo "range:" >> $output_dir/GeneralStatistics.txt
echo $range  >> $output_dir/GeneralStatistics.txt

echo " " >> $output_dir/GeneralStatistics.txt
echo "notMatchThreshold:" >> $output_dir/GeneralStatistics.txt
echo $notMatchThreshold  >> $output_dir/GeneralStatistics.txt

echo " " >> $output_dir/GeneralStatistics.txt
echo " " >> $output_dir/GeneralStatistics.txt
echo "#################################################################" >> $output_dir/GeneralStatistics.txt
echo "#################################################################" >> $output_dir/GeneralStatistics.txt
