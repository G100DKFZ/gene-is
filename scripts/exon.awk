output_dir=$1
anno9_bed=$2
echo $output_dir
cat $output_dir/$anno9_bed | awk '{ n14 = split($14, t14, ",");n15 = split($15, t15, ","); for (i = 0; ++i < n14;) {print $1,t14[i],t15[i], $10=="+" ? $9 "_Exon" i : $9 "_Exon" n14-i, $10} }' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t""0""\t"$5}' > $output_dir/exonsmyFile.bed 

