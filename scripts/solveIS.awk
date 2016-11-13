# extract the run infos from the filename. Used in LAM experiment
function getRunInfos(fields,infos) {
	infos=""
	for(i=1;i<=length(fields)-1;i++)  
		infos=infos fields[i] " ";
	return infos
}

# Print the IS in tabular or BED format (BED flag [0,1])
# BED format is used in order to generate a table that can be read via external annotation tool.

function printIS(chrPrec,IS,strand,count,flag,infos,bed) {
	if (bed)
		for (i=0;i<count;i++) 
			print chrPrec,IS,IS+1,strand
	else
		print chrPrec,IS,strand,count,flag,infos	
}

BEGIN {
	split(fileName,fields,"_");
	infos="";
	infos=getRunInfos(fields,infos);
	if (!range)
		range=10
	if (!bed)
		bed=0
	if(!countThr)
		countThr=1
	errorLam="0"
#	print "Bed mode: "bed > "/dev/tty"
#       print "Range Threshold: " range  > "/dev/tty"
#	print "Count Threshold: " countThr  > "/dev/tty"
}

{
	# check if this IS belogn to a potential IS cluster
	if ((($2-IS)*($2-IS)<=range*range)&&(chrPrec==$1)) {
	    
	    # ISTable is a symbol table <KEY=chr_pos> that is has 1 if the IS is found only 
	    # on one strand and 2 if is found on both

	    ISTable[$1"_"$2]++
            
	# countTable is a symbol table <KEY=chr_pos> that contains the sum of the reads
        # found on both strands

	    countTable[$1"_"$2]+=$4
	    
	    
	# If the method is AGILEN simply point the IS location on the IS that has the 
	# highest read count on both strands

	    if(method!="lam") {
		    if (countTable[$1"_"$2]>max) {
			    max=countTable[$1"_"$2]
			    strandPrec=$3
			    IS=$2
			    chrPrec=$1
		    }
		    count+=$4
	    }

	# If the methosd is LAM use the old definition. We should refactor the methods
	# using an ad hoc synthetic set
	    else{
		if ($3==strandPrec) {
			#if($3=="+") {
			#	if($4<count/10) 
			#		count+=$4
			#	else {
			#		if(count>=countThr)
			#			printIS(chrPrec,IS,"+",count,"0",infos,bed)
			#		IS=$2
			#		count=$4
			#	}
			#}
			#else if($3=="-") {
			if (countTable[$1"_"$2]>max) {
				max=countTable[$1"_"$2]
				strandPrec=$3
				IS=$2
				chrPrec=$1
			}
				#if($4>precCount*10) {
				#	count+=$4
				#	IS=$2
				#	precCount=$4
				#}
				#else {
				#	count+=$4
				#}
			count+=$4
		}
		else {
			#if(count>=countThr)
			printIS(chrPrec,IS,strandPrec,count,errorLam,infos,bed)
			errorLam="strand_error"
			#printIS(chrPrec,IS,strandPrec,count,"strand_error",infos,bed)
			#print chrPrec,IS,strandPrec,count,infos
			IS=$2
			count=$4
			strandPrec=$3
		}
	    }
	}
	else {
		# update the pointers and counters to the new location
		if(first==0) 
			if(count>=countThr)
				if (method!="lam") {
					if (ISTable[chrPrec"_"IS]>1)
						printIS(chrPrec,IS,strandPrec,count,"B",infos,bed)
					else printIS(chrPrec,IS,strandPrec,count,"0",infos,bed)
				}
				else printIS(chrPrec,IS,strandPrec,count,errorLam,infos,bed)
		#print chrPrec,IS,strandPrec,count,infos
		ISTable[$1"_"$2]++ 
		countTab[$1"_"$2]+=$4
		errorLam="0"
		# is the first line? (1 yes/0 no)
		first=0

		# actual chromosome 
		chrPrec=$1

		# actual location
	        IS=$2

		# actual reads count
		count=$4

		# actual strand 
		strandPrec=$3

		# max found read count
		max=$4
	}
}

END {
	if(count>=countThr)
		if (method!="lam") {
			if (ISTable[chrPrec"_"IS]>1)
				printIS(chrPrec,IS,strandPrec,count,"B",infos,bed)
			else printIS(chrPrec,IS,strandPrec,count,"0",infos,bed)
		}
		else printIS(chrPrec,IS,strandPrec,count,"0",infos,bed)

}
