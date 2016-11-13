outputDir=$1
scriptDir=$2
len=$3
alScore=$4
idClus=$5
range=$6
notMatchThreshold=$7
blatAligner=$8 
blatIndex=$9
minIden=${10}

echo "LAM sam alignment file filtering based on notMatch threshold" ;
echo " " ;
cd $outputDir
for f in *.sam
do

#		echo "for F reads"
		#get all reads with S in the start as expected and S less than 5 for flag 0 and 256
		awk  '($2=="0" || $2=="256")' $f | awk '{ split($6,a , "S");print (a[1], "\t", $0) ; }' | awk '$1  !~ /M/' | awk '$1  !~ /D/' | awk '$1  !~ /I/' | awk -v notMatchThreshold=$notMatchThreshold  '($1<=notMatchThreshold)' > ${f}_tF1
		awk  '($2=="0" || $2=="256")' $f | awk '{ split($6,a , "S");print (a[1], "\t", $0) ; }' | awk '$1  ~ /M/' > ${f}_tF2
		cat ${f}_tF1 ${f}_tF2 > ${f}_tF3
		awk '($3=="256")' ${f}_tF3 | cut -f2 | sed 's/ //g' >  ${f}_tF4
		grep -Fwf ${f}_tF4 ${f}_tF3 | cut -f2- | sed 's/ //1' | awk '($2=="0")' | cut -f1 | sed 's/ //g' > ${f}_tF5	
		grep -v -Fwf  ${f}_tF5 ${f}_tF4 > ${f}_tF6
		grep -v -Fwf   ${f}_tF6  ${f}_tF3 | cut -f2- | sed 's/ //1'  >  ${f}.tF7
		rm *_tF*

#		echo "for R reads"
		#for reads with flag 16 and 272 ...get reads with S at the end and less than 5 S part
		awk   '($2=="16" || $2=="272")'  $f   |  awk '{print substr($6,length($6),1),  "\t", $0}' | awk '($1=="S")' | cut -f2- | awk '{print substr($6,length($6)-3 ,length($6))  "\t", $0}' |  awk '{$1=A[split($1,A,"M")]}1' | sed 's/ /\t/g' |awk '{$1=A[split($1,A,"I")]}1' | sed 's/ /\t/g' | awk '{$1=A[split($1,A,"D")]}1' | sed 's/ /\t/g'  | sed 's/S//1'  | awk -v notMatchThreshold=$notMatchThreshold  '($1<=notMatchThreshold)' | cut -f2- | sed 's/ //g' |  cut -f1- > ${f}_tR1
		awk '($2=="16" || $2=="272")'  $f |  awk '{print substr($6,length($6),1),  "\t", $0}' | awk '($1=="M")' | cut -f2- | sed 's/ //g' |  cut -f1- >  ${f}_tR2
		cat ${f}_tR1 ${f}_tR2 > ${f}_tR3
		awk '($2=="272")' ${f}_tR3 | cut -f1 | sed 's/ //g' >  ${f}_tR4
		grep -Fwf ${f}_tR4 ${f}_tR3 |  awk '($2=="16")' | cut -f1  > ${f}_tR5
		grep -v -Fwf  ${f}_tR5 ${f}_tR4 > ${f}_tR6
		grep -v -Fwf   ${f}_tR6  ${f}_tR3   >  ${f}.tR7
		cat ${f}.tF7  ${f}.tR7 >  ${f}.tFR7
#		echo "rename sam file to original"
		mv  ${f} ${f}.origSam
#		echo "rename modified sam file to sam"
		mv  ${f}.tFR7  ${f}
		rm *_tR* *.tR* *.tF* 
done
#echo "................"
#echo "LAM filtering based on Threshold is finished "
#echo "................"


echo "Start IS extraction"
cd $outputDir
for f in *.sam
do
	echo "Extracting ${f}"  
        echo "awk -F\"\t\" -v ml=$len -v alScore=$alScore -f $scriptDir/lamExtractIS.awk $f"	
        awk -F"\t" -v ml=$len -v alScore=$alScore -f $scriptDir/lamExtractIS.awk $f	
	cat *.OKtoMerge  > $f"_"$len".is.repeats"
	rm *.OKtoMerge
done
echo "End IS extraction"
echo "IS Processing"
ls *.is > /dev/null 
if [ $? -eq 0 ] ;  then
	echo " Files that contain raw IS have found"
else
	echo " No files found that contains IS"
	exit
fi

#For further refining potential IS reads selected by BWA.
#Re-aligned with BLAT, and for MH compared primary and sec. hits for score and percentage identity.
#New identified MH by BLAT are added to already created MH file per sample by BWA and removed from previous potential IS file.
#Both processed further separately
cat *.trimmed2.sam > mixed
for f in *".is"
do
                IN="${f}"
                IFS='_' read -ra fields <<< "${IN}"
		cut -d " " -f1,7  $f | sed 's/^/>/'  | sed 's/ /\n/g'  > ${f}_MH1.fa
		$blatAligner  $blatIndex  ${f}_MH1.fa  -out=blast8 -minIdentity=$minIden  ${f}_MH1.bst
		sort -k1,1 -k12,12nr ${f}_MH1.bst  > ${f}_MH1.bst1a
		sed 's/\t/@/g' ${f}_MH1.bst1a > ${f}_MH1.bst1
		sort -k1,1  -u ${f}_MH1.bst1a | sed 's/\t/@/g' > ${f}_MH1.bst2
		grep -A1 -Fwf ${f}_MH1.bst2 ${f}_MH1.bst1 | sed 's/@/\t/g' | sed 's/--//g' | sed '/^$/d' > ${f}_MH1.bst3
		sort -k1,1  ${f}_MH1.bst3 | cut -f1 | uniq -d   > ${f}_MH1.bst4
		grep -Fwf ${f}_MH1.bst4 ${f}_MH1.bst3 > ${f}_MH1.bst5
		awk '{getline b;printf("%s %s\n",$0,b)}' ${f}_MH1.bst5 | sed 's/ /\t/g'| cut -f 1,2,3,4,12,15,16,24 > ${f}_MH1.bst6
		awk '{ print  $0, ($6/$3),  ($7/$4), ($8/$5)}'  ${f}_MH1.bst6 | sed 's/ /\t/g' | awk -v alScore=$alScore '($9>=alScore && $10>=alScore && $11>=alScore)' > ${f}_MH1.bst7
		cut -f1  ${f}_MH1.bst7 >  ${f}_MH1.bst8
		grep -v -Fwf ${f}_MH1.bst8  $f > ${f}_MH1.bst9
		grep  -Fwf ${f}_MH1.bst8 mixed >> ${f}.repeats
		mv  ${f}_MH1.bst9 $f
		
done
rm mixed
rm *_MH1.*



for f in *".is"
do
		IN="${f}"
		IFS='_' read -ra fields <<< "${IN}"
		lamDir=${fields[6]}
		echo "LAM direction: $lamDir"
		
		# produce IS loactions ----------------------------
		echo "cut -d \" \" -f 1,2,3,4,6,7 $f |awk -v lamDir=$lamDir 'BEGIN{if(toupper(lamDir)==\"5LAM\" || toupper(lamDir)==\"5NRLAM\") {plus=\"-\";minus=\"+\"} else {plus=\"+\";minus=\"-\"}print \"Plus: \"plus,\";Minus: \"minus > \"/dev/stderr\"}{if($4==0) $4=plus;else $4=minus;chr[$2\"\t\"$3\"\t\"$4]++}END{for (key in chr) print key\"\t\"chr[key]}'|sort -k1,1 -k2,2nr >${f}_total"  
#		cut -d " " -f 1,2,3,4,6,7 $f |awk -v lamDir=$lamDir 'BEGIN{if(toupper(lamDir)=="5LAM" || toupper(lamDir)=="5NRLAM") {plus="-";minus="+"} else {plus="+";minus="-"} print "Plus: "plus,";Minus: "minus > "/dev/stderr"}{if($4==0) $4=plus;else $4=minus;chr[$2"\t"$3"\t"$4]++}END{for (key in chr) print key"\t"chr[key]}'|sort -k1,1 -k2,2nr >${f}_total  
                cut -d " " -f 1,2,3,4,6,7 $f |awk -v lamDir=$lamDir 'BEGIN{if(toupper(lamDir)=="5LAM" || toupper(lamDir)=="5NRLAM") {plus="-";minus="+"} else {plus="+";minus="-"}}{if($4==0) $4=plus;else $4=minus;chr[$2"\t"$3"\t"$4]++}END{for (key in chr) print key"\t"chr[key]}'|sort -k1,1 -k2,2nr >${f}_total
	
		echo "cat ${f}_total | sort -k1,1 -k2,2n | awk  -v method=lam -v fileName=$f -v range=$range -f $scriptDir/solveIS.awk > ${f}.IS.results"
		cat ${f}_total | sort -k1,1 -k2,2n | awk -v method=lam -v fileName=$f -v range=$range -f $scriptDir/solveIS.awk > ${f}.IS.results
		
		# Produce clustered sequences ------------------------
		echo "cut -d\" \" -f 1,2,3,4,6,7 $f|awk 'BEGIN{plus="+";minus="-"}{if(toupper(lamDir)=="5LAM") {plus="-";minus="+"} if ($4==0) $4=plus; else $4=minus;print $0}' |awk '{printf(\">\");for (i=2;i<=NF;i++) printf(\"%s \",$i);printf(\"%s \",$1);print \"\n\"$NF}'> ${f}_fa"  
		cut -d" " -f 1,2,3,4,6,7 $f|awk 'BEGIN{plus="+";minus="-"}{if(toupper(lamDir)=="5LAM") {plus="-";minus="+"}if ($4==0) $4="+"; else $4="-";print $0}' |awk '{printf(">");for (i=2;i<=NF;i++) printf("%s ",$i);printf("%s ",$1);print "\n"$NF}'> ${f}_fa
		
###		echo "$scriptDir/uclust --sort ${f}_fa --output ${f}_sort"
###		$scriptDir/uclust --sort ${f}_fa --output ${f}_sort   
###		echo "$scriptDir/uclust --input ${f}_sort --rev --id $idClus --uc ${f}_res"
###		$scriptDir/uclust --input ${f}_sort --rev --id $idClus --uc ${f}_res  
		
###		awk '$1=="C"' ${f}_res |sort -k3nr|cut -f 3,9,10,11,12,13 > ${f}.collapsed
###		echo "point6"
		# Merging IS table with collapsed sequences -------------------
		#join -a 1 <(awk '{print $1"_"$2,$0}' ${f}.IS.results|sort -k1,1) <(sort -k2,2 -k3,3n ${f}.collapsed|awk '{if(pos==$3) {if ($1>max) {max=$1; out=$6}} else {print chr"_"pos,max,out;chr=$2;pos=$3;max=$1;out=$6}}END{print chr"_"pos,max,out}'|sort -k1,1))|awk '{$1="";print$0}' > ${f}.seq
###		join -a 1 <(awk '{print $1"_"$2,$0}' ${f}.IS.results|sort -k1,1) <(sort -k2,2 -k3,3n ${f}.collapsed|awk -f  $scriptDir/getClusterSequences.awk|sort -k1,1)|awk '{$1="";print$0}' > ${f}.seq	
###		echo "cat MiS*.seq > resultsNoDup.csv.total.results"
###		cat MiS*.seq > resultsNoDup.csv.total.results	

		echo "cat MiS*.IS.results > resultsNoDup.csv.total.results"
                cat MiS*.IS.results > resultsNoDup.csv.total.results

done
