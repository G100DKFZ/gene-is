id=0;
tagTable=$1
output_dir=$2
script_dir=$3
threads=$4
readHeader=$5

cd $output_dir/tmpR;
for f in *;
do  
	echo "Started Sorting $id";
	echo "awk -v fileIn=$tagTable -v outDir=$output_dir -v procID=$id -v trimming=0 -f $script_dir/SortLtrLnkParallelTest.awk $f $output_dir/tmpR/$f &"
	awk -v fileIn=$tagTable -v outDir=$output_dir -v procID=$id -v trimming=0 -v readHeader=$readHeader -f $script_dir/SortLtrLnkParallelTest.awk $output_dir/tmpF/$f $output_dir/tmpR/$f &
	(( ++id%threads == 0 )) && wait;
done;
wait
echo "cd ${output_dir}"
cd ${output_dir}
#echo "ls | awk -F\"__\" '{cmd=\"cat \" $0 \" >> \" $1 \"_.fastq\"; system (cmd)}'"
ls | awk -F"__" '{cmd="cat " $0 " >> " $1 "_.fastq"; system (cmd)}'
echo "rm *__*;rm -r tmp*"
rm *__*;rm -r tmp*
