output_dir=$1
indexFile=$2
aligner=$3
threads=$4
samtools=$5
echo "Start Alignment Phase"
echo "on $output_dir"
cd $output_dir
for f in *trimmed2
do
	IN="${f}"
	IFS='_' read -ra fields <<< "${IN}"
	lamDir=${fields[6]}
	echo "Aligning $f" >&2
	echo "LAM direction: $lamDir" >&2
	if [[ $aligner == *bwa* ]]; then
		echo "$aligner mem -t $threads -M $indexFile  $f  > ${f}.sam" >&2
		$aligner mem -M -t $threads $indexFile  $f   > ${f}.sam
		if [[ $lamDir == *nrLAM* ]]; then
			echo "mv ${f}.sam swap.sam" >&2
			mv ${f}.sam swap.sam	
			echo "$samtools view -Sb swap.sam > swap.bam" >&2
			$samtools view -Sb swap.sam > swap.bam
	        	echo "$samtools sort swap.bam swap.sorted" >&2
	        	$samtools sort swap.bam swap.sorted
			echo "$samtools index swap.sorted.bam" >&2
			$samtools index swap.sorted.bam
			echo "$samtools rmdup -s swap.sorted.bam swap.nodup.bam" >&2
			$samtools rmdup -s swap.sorted.bam swap.nodup.bam
			#echo "mv ${f}.sam ${f}.sam.dup" >&2
			#mv ${f}.sam ${f}.sam.dup
			echo "$samtools view -h -o swap.sam swap.nodup.bam" >&2
			$samtools view -h -o swap.sam swap.nodup.bam 
			echo "mv swap.sam ${f}.sam" >&2
			mv swap.sam ${f}.sam
		fi
	else
		echo "$aligner -p $threads --very-sensitive-local -x $indexFile -U $f -S ${f}.sam"
		$aligner -p $threads --very-sensitive-local -x $indexFile -U $f -S ${f}.sam 2>/dev/null
		if [[ $lamDir == *nrLAM* ]]; then 
			echo "$samtools view -Sb ${f}.sam > ${f}.bam"
			$samtools view -Sb ${f}.sam > ${f}.bam
	        	echo "$samtools sort ${f}.bam ${f}.sorted"
	        	$samtools sort ${f}.bam ${f}.sorted
			echo "$samtools index ${f}.sorted.bam"
			$samtools index ${f}.sorted.bam
			echo "$samtools rmdup -s ${f}.sorted.bam ${f}.nodup.bam"
			$samtools rmdup -s ${f}.sorted.bam ${f}.nodup.bam
			echo "mv ${f}.sam ${f}.sam.dup"
			mv ${f}.sam ${f}.sam.dup
			echo "$samtools view -h -o ${f}.sam ${f}.nodup.bam"
			$samtools view -h -o ${f}.sam ${f}.nodup.bam
		fi
	fi
done
echo "End Alinment Phase"
