outputDir=$1
scriptDir=$2
len=$3
alScore=$4
idClus=$5
range=$6

echo "Start Multiple Aligned IS extraction"
cd $outputDir

rename 's/.is.repeats$/.is..r/' *.is.repeats
less *.is |  cut -d " " -f1  > tt1
less *.is..r |  cut  -f1  > tt2
grep -Fwf tt2 tt1 > tt3

for f in *".is..r"
do
                IN="${f}"
                IFS='_' read -ra fields <<< "${IN}"
                grep -v -Fwf  tt3  $f  >  ${f}.is.r
done



for f in *.is.r
do
	#echo "Extracting repeats/multipleAligned  ${f}"  
        echo "awk -F\"\t\" -v ml=$len -v alScore=$alScore -f $scriptDir/repeatsExtractLAM.lamExtractIS.awk $f"	
        awk -F"\t" -v ml=$len -v alScore=$alScore -f $scriptDir/repeatsExtractLAM.lamExtractIS.awk  $f	
done
echo "End  Multiple Aligned IS extraction"
echo " Multiple Aligned IS Processing"
ls *.is > /dev/null 
if [ $? -eq 0 ] ;  then
	echo " Files that contain raw IS repeats/multipleAligned have found"
else
	echo " No files found that contains IS repeats/multipleAligned"
	exit
fi
for f in *".is.r1"
do
		IN="${f}"
		IFS='_' read -ra fields <<< "${IN}"
		lamDir=${fields[6]}
		echo "LAM direction: $lamDir"

		# produce IS loactions ----------------------------
		echo "cut -d \" \" -f 1,2,3,4,6,7 $f |awk -v lamDir=$lamDir 'BEGIN{if(toupper(lamDir)==\"5LAM\" || toupper(lamDir)==\"5NRLAM\") {plus=\"-\";minus=\"+\"} else {plus=\"+\";minus=\"-\"}print \"Plus: \"plus,\";Minus: \"minus > \"/dev/stderr\"}{if($4==0) $4=plus;else $4=minus;chr[$2\"\t\"$3\"\t\"$4]++}END{for (key in chr) print key\"\t\"chr[key]}'|sort -k1,1 -k2,2nr >${f}_total.r"  
#		cut -d " " -f 1,2,3,4,6,7 $f |awk -v lamDir=$lamDir 'BEGIN{if(toupper(lamDir)=="5LAM" || toupper(lamDir)=="5NRLAM") {plus="-";minus="+"} else {plus="+";minus="-"} print "Plus: "plus,";Minus: "minus > "/dev/stderr"}{if($4==0) $4=plus;else $4=minus;chr[$2"\t"$3"\t"$4]++}END{for (key in chr) print key"\t"chr[key]}'|sort -k1,1 -k2,2nr >${f}_total.r
	cut -d " " -f 1,2,3,4,6,7 $f |awk -v lamDir=$lamDir 'BEGIN{if(toupper(lamDir)=="5LAM" || toupper(lamDir)=="5NRLAM") {plus="-";minus="+"} else {plus="+";minus="-"}}{if($4==0) $4=plus;else $4=minus;chr[$2"\t"$3"\t"$4]++}END{for (key in chr) print key"\t"chr[key]}'|sort -k1,1 -k2,2nr >${f}_total.r	        
		echo "cat ${f}_total.r | sort -k1,1 -k2,2n | awk  -v method=lam -v fileName=$f -v range=$range -f $scriptDir/solveIS.awk > ${f}.IS.results.r"
		cat ${f}_total.r | sort -k1,1 -k2,2n | awk -v method=lam -v fileName=$f -v range=$range -f $scriptDir/solveIS.awk > ${f}.IS.results.r
		
		# Produce clustered sequences ------------------------
		echo "cut -d\" \" -f 1,2,3,4,6,7 $f|awk 'BEGIN{plus="+";minus="-"}{if(toupper(lamDir)=="5LAM") {plus="-";minus="+"} if ($4==0) $4=plus; else $4=minus;print $0}' |awk '{printf(\">\");for (i=2;i<=NF;i++) printf(\"%s \",$i);printf(\"%s \",$1);print \"\n\"$NF}'> ${f}_fa.r"  
		cut -d" " -f 1,2,3,4,6,7 $f|awk 'BEGIN{plus="+";minus="-"}{if(toupper(lamDir)=="5LAM") {plus="-";minus="+"}if ($4==0) $4="+"; else $4="-";print $0}' |awk '{printf(">");for (i=2;i<=NF;i++) printf("%s ",$i);printf("%s ",$1);print "\n"$NF}'> ${f}_fa.r
		
###		echo "$scriptDir/uclust --sort ${f}_fa.r --output ${f}_sort.r"
###		$scriptDir/uclust --sort ${f}_fa.r --output ${f}_sort.r 
###		echo "$scriptDir/uclust --input ${f}_sort.r --rev --id $idClus --uc ${f}_res.r"
###		$scriptDir/uclust --input ${f}_sort.r --rev --id $idClus  --uc ${f}_res.r  
		
###		awk '$1=="C"' ${f}_res.r |sort -k3nr|cut -f 3,9,10,11,12,13 > ${f}.collapsed.r
		
		# Merging IS table with collapsed sequences -------------------
		#join -a 1 <(awk '{print $1"_"$2,$0}' ${f}.IS.results.r |sort -k1,1) <(sort -k2,2 -k3,3n ${f}.collapsed.r |awk '{if(pos==$3) {if ($1>max) {max=$1; out=$6}} else {print chr"_"pos,max,out;chr=$2;pos=$3;max=$1;out=$6}}END{print chr"_"pos,max,out}'|sort -k1,1))|awk '{$1="";print$0}' > ${f}.seq.r
###		join -a 1 <(awk '{print $1"_"$2,$0}' ${f}.IS.results.r |sort -k1,1) <(sort -k2,2 -k3,3n ${f}.collapsed.r |awk -f  $scriptDir/getClusterSequences.awk|sort -k1,1)|awk '{$1="";print$0}' > ${f}.seq.r	
###		echo "cat MiS*.seq.r > resultsNoDup.csv.total.results.r"
###		cat MiS*.seq.r  > resultsNoDup.csv.total.results.r	

                echo "cat MiS*.IS.results.r > resultsNoDup.csv.total.results.r"
                cat MiS*.IS.results.r  > resultsNoDup.csv.total.results.r

done
