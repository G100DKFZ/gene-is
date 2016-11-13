# This script takes all the fastq files in the output directory (first argument provided) and remove all the sequencesthat contain:
# - the megaprimer (present in the column 18 of the tagTable [field 17 of the fileanem])
#   If the megaprimer is not present (len(sequence)!=len(quality)) the read is discarded.
#   Only perfect sequences are taken into account.
# - the linker (present in the column 21 of the tagTable [field 20 of the fileanem])
# - The restirction enzyme id taken from field 8 (7 in filename)
# - a sequence added from the sequencer that is flagged by the pattern [TTGAAAAAA.*] (should be better modelled)
# 

# Hash table for restriction enzymes; please fill with all the enzymes in use
# for enzymes in witch the restriction site is changed (due to the addiction of the barcode, ex MseI (TTAA -> TTATAG) fake restriction sites are introduced...be aware)
# declare -A enzymes=( ["MseI"]="TTAA" ["FatI"]="CATG" ["HpyCH4IV"]="ACGT" ["MluCI"]="AATT" ["Tsp509I"]="AATT" ["circLigase"]="circligase")

#megaprimer=$1
#linker=$2
output_dir=$1
suff=$2
qual=$3
skewer_dir=$4
#echo "cd ${output_dir}"
cd ${output_dir}
#ls | awk -F"__" '{cmd="cat " $0 " >> " $1 ".fastq"; system (cmd)}'
for f in $suff*.fastq;
do 
	echo \"Trimming ${f}\";
	IN="${f}"
	IFS='_' read -ra fields <<< "${IN}"
	megaprimer=${fields[17]}
	linker=${fields[20]}
	enzyme=${fields[7]}
	bcLinker=${fields[22]}
	echo "Megaprimer-> $megaprimer"
	echo "Restriction Enzyme: $enzyme -> ${enzymes[${enzyme}]}"
	
	
	# uncomment next line if you want to remove linker string and poly A stretch introduced by the sequencer when len(amplicon)<len(read)
	# sed -e "s/^.*${megaprimer}//" $f | awk -v outputDir=$output_dir '{id=$0; getline; seq=$0; len=length($0);getline;getline;if (len!=length($0)) print id"\n"seq"\n+\n"substr($0,length($0)-len+1); else print id "\n" seq "\n+\n" $0 >> outputDir"/LeftNoMega.txt"}'| sed -e  "s/${enzymes[${enzyme}]}.*/${enzymes[${enzyme}]}/" |sed -e  "s/${linker}.*//" |sed -e  "s/T[TG][CG]AAAAAA.*//" | awk '{print $0; getline; print $0; len=length($0);getline;print $0;getline;print substr($0,1,len)}' > ${f}.trimmed2
	
	
	# comment following line in order to do not use restriction enzyme detection
	#sed -e "s/^.*${megaprimer}//" $f | awk -v outputDir=$output_dir '{id=$0; getline; seq=$0; len=length($0);getline;getline;if (len!=length($0)) print id"\n"seq"\n+\n"substr($0,length($0)-len+1); else print id "\n" seq "\n+\n" $0 >> outputDir"/LeftNoMega.txt"}'| sed -e  "s/${enzymes[${enzyme}]}.*/${enzymes[${enzyme}]}/" | awk '{print $0; getline; print $0; len=length($0);getline;print $0;getline;print substr($0,1,len)}' > ${f}.trimmed2

	# Use skewer to trimm linker sequence and bad quality (<30) read portion	
	echo "$skewer_dir -x GGGGGGGGGGGGGGGGGGGG -d 0 -q $qual -m tail -o test ${f}"
	$skewer_dir -x GGGGGGGGGGGGGGGGGGGG -d 0 -q $qual -m tail -o test ${f}
	echo "$skewer_dir -x $bcLinker -l 25 -b -r 0.2 -m tail -o ${f} test.fastq"
	$skewer_dir -x $bcLinker -l 25 -b -r 0.2 -m tail -o ${f} test.fastq
	mv ${f}X01.fastq ${f}.link
	echo $#;
	# Megaprimer portion (sequence from tag table) is detected with 100% identity. If not detected the read is discarded 
	if [ $# -eq 4 ]; then
		sed -e "s/^.*${megaprimer}//" $f.link | awk -v outputDir=$output_dir '{id=$0; getline; seq=$0; len=length($0);getline;getline;if (len!=length($0)) print id"\n"seq"\n+\n"substr($0,length($0)-len+1); else print id "\n" seq "\n+\n" $0 >> outputDir"/LeftNoMega.txt"}' > ${f}.trimmed2
	fi
	# An anchor is given to the trimming process
	if [ $# -eq 6 ] ; then
		anchor=$5
		anchorMM=$6
		sed -e "s/^.*${megaprimer}//" $f.link | awk -v outputDir=$output_dir '{id=$0; getline; seq=$0; len=length($0);getline;getline;if (len!=length($0)) print id"\n"seq"\n+\n"substr($0,length($0)-len+1); else print id "\n" seq "\n+\n" $0 >> outputDir"/LeftNoMega.txt"}' > ${f}.mega
		echo "$skewer_dir -x $anchor -d $anchorMM -l 25 -b -m head -o anchor ${f}.mega"
		$skewer_dir -x $anchor -d $anchorMM -l 25 -b -m head -o anchor ${f}.mega
		echo "mv anchorX01.fastq ${f}.trimmed2"
		mv anchorX01.fastq ${f}.trimmed2
		#sed -e "s/^.*${megaprimer}//" $f.link | awk -v outputDir=$output_dir '{id=$0; getline; seq=$0; len=length($0);getline;getline;if (len!=length($0)) print seq"\t"substr($0,length($0)-len+1)"\t"id; else print id "\n" seq "\n+\n" $0 >> outputDir"/LeftNoMega.txt"}' > ${f}.mega
		#agrep -$anchorMM "^$anchor" ${f}.mega | awk -F"\t" '{print $3"\n"$1"\n+\n"$2}' > ${f}.trimmed2
	fi
done

