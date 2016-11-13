output_dir=$1
anno12_bed=$2
ISFileMod0_bed=$3

cut -f1,2,3,4,5,6,7,8 $output_dir/$anno12_bed  > $output_dir/anno13.bed
cut -f 9,10,11,12 $output_dir/$anno12_bed | sed 's/-/ /g'  | sed 's/+/ /g' > $output_dir/anno14.bed
paste $output_dir/anno13.bed $output_dir/anno14.bed > $output_dir/anno15.bed
echo -e  "Chr\tIS\tStrand\tSeq_Count\tRefSeq_ID\tGene_Strand\tGene_Name\tGene_Length\tDist_to_TSS\tUpstream\tDownstream\tIntron_Exon" > $output_dir/result_AnnotatedIS.csv 
awk '($8!="0")'  $output_dir/anno15.bed > $output_dir/anno16.bed
awk '($8=="0")' $output_dir/anno15.bed  | awk '{print ($1"\t"$2"\t"$3"\t"$4"\t""NA""\t""NA""\t""NA""\t""NA""\t""NA""\t""NA""\t""NA""\t""NA""\t")}' >  $output_dir/anno17.bed
cat $output_dir/anno16.bed $output_dir/anno17.bed > $output_dir/anno18.bed
#addition to cut all columns from LAM and cat to output
sort -k1,1 -k2,2n  -k3,3 -k4,4n   $output_dir/$ISFileMod0_bed >  $output_dir/ISFileMod00.bed
sort -k1,1 -k2,2n  -k3,3 -k4,4n  $output_dir/anno18.bed > $output_dir/anno18a.bed
cut -f7- $output_dir/ISFileMod00.bed > $output_dir/ISFileMod000.bed
#fix column shift
sed 's/\t/_/g'  $output_dir/ISFileMod000.bed >  $output_dir/ISFileMod000.0.bed
cut -f1,2,3,4,5,6,7,8,9,10,11 $output_dir/anno18a.bed > $output_dir/anno18a.a1.bed
 cut -f12 $output_dir/anno18a.bed > $output_dir/anno18a.a2.bed
 paste $output_dir/anno18a.a1.bed  $output_dir/anno18a.a2.bed  $output_dir/ISFileMod000.0.bed  | sed 's/\t/,/g' | sort  -t',' -k1,1 -k2,2 -k13,13  -u | sed 's/,/\t/g' >> $output_dir/result_AnnotatedIS.csv
sed 's/\t/,/g' $output_dir/result_AnnotatedIS.csv >  $output_dir/result_AnnotatedIS.csv1
mv $output_dir/result_AnnotatedIS.csv1 $output_dir/ResultsClusteredAnnotated.csv
#result file with ten strongest clone
echo -e  "Chr\tIS\tStrand\tSeq_Count\tRefSeq_ID\tGene_Strand\tGene_Name\tGene_Length\tDist_to_TSS\tUpstream\tDownstream\tIntron_Exon" >  $output_dir/ResultsTenStrongestClones.csv
tail -n +2  $output_dir/result_AnnotatedIS.csv  |  sort -rn -k4,4  | head -10  >> $output_dir/ResultsTenStrongestClones.csv
sed -i 's/\t/,/g' $output_dir/ResultsTenStrongestClones.csv
rm $output_dir/*bed
