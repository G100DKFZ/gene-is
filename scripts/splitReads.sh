id=0;
forward=$1
reverse=$2
tagTable=$3
suff=$4
lenFilter=$5
threads=$5
group=$6
megaprimer=$7
linker=$8
AssemblyDir=$9
BowtieIndx=$10


# This part should have the signature
# -f <fastq.gz>
# -r <fastq.gz>
# -g <int> [5000000]
echo "Forward Splitting started"

time zcat $forward| split -l $group - tmpR/test &



######time zcat $1| head -1000000| split -l 100000 - tmpR/test &
######zcat $1| head -12| split -l 1000 - tmpR/test &



echo "Reverse Splitting started"
time zcat $reverse|split -l $group - tmpL/test &


######time zcat $2|head -1000000| split -l 100000 - tmpL/test &
######zcat $2|head -12| split -l 1000 - tmpL/test &
wait
echo "Stop Splitting"
cd tmpR
for f in test*;
do
	echo "Started Sorting $id"
	#time awk -v fileIn=$tagTable -v suff=$suff -v procID=$id -f ../SortLtrLnkParallelTest.awk $f ../tmpL/$f &
	#(( ++id%threads == 0 )) && wait;
done
echo "End Sorting\nStarted Merging"
cd ../Sorted${suff}
#ls | awk -F"__" '{cmd="cat " $0 " >> " $1 ".fastq"; system (cmd)}'

echo "End Merging"
echo "Start Trimming"
for f in MiS43_MiSeqPaper*fastq 
do
        echo "Trimming ${f}"	
	#skewer -x $megaprimer -l 25 -m head -q 30 -o ${f}trimmed $f 
	#skewer -x GATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA -l 25 -m head -q 30 -o ${f}trimmed $f 
	#sed -e 's/^.*${megaprimer}//' $f |awk '{id=$0; getline; seq=$0; len=length($0);getline;getline;if (len!=length($0)) print id"\n"seq"\n+\n"substr($0,length($0)-len+1)}'|sed -e  "s/CCTAACTGCTGTGCCACT.*//" ${f}trimmed.fastq|sed -e  "s/TTGAAAAAA.*//" |awk '{print $0; getline; print $0; len=length($0);getline;print $0;getline;print substr($0,1,len)}' > ${f}trimmed2



#	sed -e 's/^.*GATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA//' $f | awk '{id=$0; getline; seq=$0; len=length($0);getline;getline;if (len!=length($0)) print id"\n"seq"\n+\n"substr($0,length($0)-len+1)}'| sed -e  "s/CCTAACTGCTGTGCCACT.*//" |sed -e  "s/TTGAAAAAA.*//" | awk '{print $0; getline; print $0; len=length($0);getline;print $0;getline;print substr($0,1,len)}' > ${f}trimmed2
	
	#sed -e 's/^.*GATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA//' $f | awk '{id=$0; getline; seq=$0; len=length($0);getline;getline;if (len!=length($0)) print id"\n"seq"\n+\n"substr($0,length($0)-len+1); else print id"\n"seq"\n+\n"$0 >> "LeftNoMega.txt"}'| sed -e  "s/CCTAACTGCTGTGCCACT.*//" |sed -e  "s/TTGAAAAAA.*//" | awk '{print $0; getline; print $0; len=length($0);getline;print $0;getline;print substr($0,1,len)}' > ${f}trimmed2
done
echo "End Trimming"
echo "Start Alignment"
for f in *fastqtrimmed2* 
do 
        echo "Aligning ${f}"	
	#~/bin/bowtie2-2.2.2/bowtie2 -p 8 --very-sensitive-local -x ${AssemblyDir}/${BowtieIndx} -U $f  -k 2 -S ${f}.sam 2>/dev/null 
	#~/bin/bowtie2-2.2.2/bowtie2 -p 8 --very-sensitive-local -x ../../HISAP2/VirusFinder2.0/TestSamples/hg19+virus -U $f  -k 2 -S ${f}.sam 2>/dev/null 
done
echo "End Alignemt"
echo "Start IS extraction"
for f in *11u*fastqtrimmed2.sam 
do 
        echo "Extracting ${f}"	
	#awk -F"\t" -v ml=$lenFilter -f extractIS.awk $f
	awk -F"\t" -v ml=25 -f ../extractIS.awk $f

done
echo "End IS extraction"
echo "IS Processing"
for f in *11u*.is; do
	cut -d " " -f 1,2,3,4,6,7 $f |awk '{if($4==0) $4="+";else $4="-";chr[$2"\t"$3"\t"$4]++}END{for (key in chr) print key"\t"chr[key]}'|sort -k4nr >${f}_total
done
for f in *_25.is 
do 
	cut -d" " -f 1,2,3,4,6,7 $f|awk '{if($4==0) $4="+";else $4="-";print $0}' | awk '{printf(">");for (i=2;i<=NF;i++) printf("%s\t",$i);printf("%s\t",$1);print "\n"$NF}' > ${f}_fa
	uclustq1.2.22_i86linux64 --sort ${f}_fa --output ${f}_sort
	uclustq1.2.22_i86linux64 --input ${f}_sort --uc ${f}_res --id 0.95
	awk '$1=="C"' ${f}_res|sort -k3nr|cut -f 3,9,10,11,12,13 > ${f}.collapsed	
done
echo "End IS processing"

