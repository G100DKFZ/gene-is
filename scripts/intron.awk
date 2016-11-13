output_dir=$1
anno9_bed=$2

cat $output_dir/$anno9_bed | awk '{ n14 = split($14, t14, ",");n15 = split($15, t15, ",");for (i = 0; ++i < n14-1;) { m14= n14-1; print $1,t15[i],t14[i+1], $10=="+" ? $9 "_Intron" i : $9 "_Intron"  m14-i, $10 } }'  | awk '{print $1"\t"$2"\t"$3"\t"$4"\t""0""\t"$5}' >  $output_dir/intronsmyFile.bed 
