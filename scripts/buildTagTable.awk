# This awk script takes a csv TagTable of 21 columns and return a tabular file
# of 6 columns that is used in order to perform a double sorting
# 
#
#
#
#
# Usage awk -f buildTagTable.awk <TagTable>





BEGIN{
        tab["A"]="T"
        tab["C"]="G"
        tab["G"]="C"
        tab["T"]="A"
	code=0
	lnkCount=1
	ltCount=1
}{
	bc_="";ltBC="";name="";bcRev_=""
	linkerSeqRev="";linkerSeq="";tagSeq="";tagSeqRev=""
	name=$11"_"$1
	# seqrch in field 6 if is a PCR or not
	# If yes then use field 15 as templlate for reverse primer and do not use the classical Linker sequence
	if ($6=="PCR") {
		split(toupper($15),invertPCR,"")
		#for(l=length(invertPCR);l>0;l--)
		#	invertPrimer=invertPrimer tab[invertPCR[l]]

		linkerSeqRev=substr(toupper($15),1,6)
		split(linkerSeqRev,invertLink,"")
		for(l=length(invertLink);l>0;l--)
			linkerSeq=linkerSeq tab[invertLink[l]]
		print linkerSeq " " linkerSeqRev
	} else { 
		#linkerSeq="CCTAAC"
		#linkerSeqRev="GTTAGG"
		linkerSeq=substr(toupper($21),1,6)
		split(linkerSeq,invertLink,"")
		for(l=length(invertLink);l>0;l--)
			linkerSeqRev = linkerSeqRev tab[invertLink[l]]
		tagSeq="TAG"
		tagSeqRev="CTA"
	}	
        for (i=2;i<12;i++)
		name=name"_"$i
	# Complement the lkBC uncomment next 3 lines
	split($12,invertBC,"")
        for(l=length(invertBC);l>0;l--)
           bc_=bc_ tab[invertBC[l]]
	bcRev_=toupper($12)
		
	if ($6=="PCR") {
	       #print "bcrev->" bcRev_ " lkRev->" linkerSeqRev " bc->" bc_ " lk->"linkerSeq_	
		bcRev_= bcRev_ linkerSeqRev
		bc_= linkerSeq bc_
	}
	else	{
		bc_=tagSeq bc_ linkerSeq
		bcRev_=linkerSeqRev bcRev_ tagSeqRev
	}
	
	
	#if (length(bcRev)==12)
	#	bcRev="TAG" bcRev "CCTAAC"
	#else
	#	bcRev="TAG" bcRev "CCTAACTG"

	#No complement lkBC uncomment nex line
	

	#if (length(bc)==12)
	#	bc="GTTAGG" bc "CTA"
	#else
	#	bc="CAGTTAGG" bc "CTA"

	for (i=13;i<22;i++)
        	name=name"_"$i
	ltBC=toupper($17) toupper(substr($18,1,6))
	name=name "_"$12 "_" bc_
	if(ltrTable[ltBC]=="") {
		ltrTable[ltBC]="ltBC" ltCount
		ltCount++
	}
	if(lnkTable[bcRev_]=="") {
		lnkTable[bcRev_]="lkBC" lnkCount
		lnkCount++
	}
	tmpLtrBC=ltrTable[ltBC]
	tmpLnkBC=lnkTable[bcRev_]
	print ltBC "\t" bcRev_ "\t" bc_ "\t" name "\tExp" code "\t"tmpLtrBC"\t"tmpLnkBC > "tagTab"$11".txt"
	
	# Uncomment if you want two separate barcode files (1 for ltr and 1 for lnkr)
	#print "Exp"code"\t"ltBC > "ltrBCset"$11".txt"
	#print "Exp"code"\t"bc > "lnkBCset"$11".txt"
	code++
}
END {
#	for (key in ltrTable) 
#		print ltrTable[key]"\t"key>"ltrBCset"$11".txt"
#	for (key in lnkTable) 
#		print lnkTable[key]"\t"key>"lnkBCset"$11".txt"
	print "#LtrBC-> " ltCount-1 "\t#LnkBC-> " lnkCount
}	
