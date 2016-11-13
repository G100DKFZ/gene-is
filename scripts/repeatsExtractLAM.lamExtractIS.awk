function extractMatch (stringa,flag,idx) {
	split(stringa,supp,flag)
	split(supp[idx],supp2,"M")
	return supp2[1]
}

BEGIN{
	#LTRfragLen=20
	#LTRfragLen=37
	LTRfragLen=155
}
{
	# alScore define the maximal homology ratio ]0,1] allowed between primary and secondary BWA aligment
	# see Calabria a et al, 2014 for details
####	if (alScore>0) {
		# AS:i:value is the BWA field that ocntains the score of the primary alignment
####		if(pos=match($0,"AS:i:")) {
####			split(substr($0,pos),vect,"\t")
####			split(vect[1],AS,":")
			# check if the read is aligned
####			if (AS[3]>0) {
				# XS:i:val is the BWA field that contains the score of the secondary alignment (0 if not present)
####				split(vect[2],XS,":")
####				if (XS[3]/AS[3] >= alScore) {
####					print $0 > fileName[1]"_"ml".is.repeats"
####					#discard the alignment
####					next
####				}
####			}
####		}
####	}
	if ((match($3,"chr"))&&($2<=16)&&(length($10)>=ml)) {
		split(FILENAME,fileName,"fastq")
		split($6,cigar1,"S")
		# If the softclip is not present the read is malformed. Uncomment the nex two lines
		# to remove this filtering step
		#if (length(cigar1)==1) {
		#	print $0 "-> Complete Vector Deletion"  > fileName[1]"_"ml".is.removed"
		#	next;
		#}
		
		
		split($6,cigar2,"M")
		# Caso CIGAR=XMYIZM
		# The read is Malformed because a deletion of LTRfragLen bases (in caseyy of typical megaprimer) produce this kind of case
		if ((length(cigar1)==1)&&(length(cigar2)==3)&&((match($6,"I"))||(match($6,"D")))) {
			# comment the next 2 lines if you want also to consider this case
			#print $0 "-> Complete Vector Deletion"  > fileName[1]"_"ml".is.removed"
			#next;
			
			
			if (match($6,"I"))
			       flag="I"
			else if (match($6,"D"))
				flag="D"
			numMatch=extractMatch(cigar1[1],flag,1)+extractMatch(cigar1[1],flag,2);
			if($2==16) {
				delta=numMatch
			}
			else {
				delta=0;
			}
			if ( numMatch>=ml)

				print $1" "$3" "$4+delta" "$2" "$6" " numMatch " "$10 > fileName[1]"_"ml".is.r1"
				#print $1" "$3" "$4+delta" "$2" "$6" " numMatch " "$10 
			else 
				 print $0 "-> Seq length smaller than " ml   > fileName[1]"_"ml".is.removed.r"
			next;
		}
		# Caso 1: CIGAR=XMYIZMKS
		# Caso 2: CIGAR=kSXMYIZM
		if ((length(cigar1)==2)&&(length(cigar2)==3)&&((match($6,"I"))||(match($6,"D")))) {
			if (match($6,"I"))
			       flag="I"
			else if (match($6,"D"))
				flag="D"
			cigar1[2]==""?caso=1:caso=2;
			
			# If the 36 (+-4) bases are not found the IS is removed
		        #if (cigar1[caso]<LTRfragLen-150 ||cigar1[caso]>LTRfragLen+4) {
			#	print $0 "-> Deletion smallest than admitted ("cigar1[caso]" < " LTRfragLen-150")"  > fileName[1]"_"ml".is.removed"
			#	next;
			#}
			
			
			numMatch=extractMatch(cigar1[caso],flag,1)+extractMatch(cigar1[caso],flag,2);
			if($2==16) {
				delta=numMatch
				sequenza=substr($10,1,length($10)-cigar1[caso])
			}
			else {
				delta=0;
				sequenza=substr($10,cigar1[caso])
			}
			
			if ( numMatch>=ml)
				print $1" "$3" "$4+delta" "$2" "$6" " numMatch " "sequenza > fileName[1]"_"ml".is.r1"
				#print $1" "$3" "$4+delta" "$2" "$6" " numMatch " "sequenza
			else 
				print $0 "-> Seq length smaller than " ml   > fileName[1]"_"ml".is.removed.r"
			next;
		}

		# Caso : CIGAR=XSYMXS
		if(length(cigar1)==3 && length(cigar2)==2) {
			split(cigar1[2],IS1,"M")
			$2==16?delta=IS1[1]:delta=0
			#print IS1[1]
			if(IS1[1]>=ml)
				print $1" "$3" "$4+delta" "$2" "$6" "IS1[1]" "$10 > fileName[1]"_"ml".is.r1"
			else 
				print $0 "-> Seq length smaller than " ml   > fileName[1]"_"ml".is.removed.r"
			next
		}

		# Caso : CIGAR=XM
		# The read is Malformed because a deletion of 22 bases (in caseyy of typical megaprimer) produce this kind of case
		if (length(cigar1)==1) {
				# comment the next 2 lines if you want also to consider this case
				#print $0 "-> Complete Vector Deletion"  > fileName[1]"_"ml".is.removed"
				#next;
				$2==16?delta=cigar2[1]:delta=0
				if(cigar2[1]>=ml)
					print $1" "$3" "$4+delta" "$2" "$6" "cigar2[1]" "$10 > fileName[1]"_"ml".is.r1"
					#print $1" "$3" "$4+delta" "$2" "$6" "cigar2[1]" "$10 
				else 
					print $0 "-> Seq length smaller than " ml   > fileName[1]"_"ml".is.removed.r"
				next
		}
		# Caso : CIGAR=XMYS or XSYM
		if ((length(cigar1)<=2)&&(length(cigar2)<=2)) {
				#print cigar1[1],cigar2[1]
				split(cigar1[1],IS1,"M")
				
				if (length(IS1)==2) {
					# XMYS
					# If the 22 (+-4) bases are not found the IS is removed
		        		#if (IS1[2]<LTRfragLen-150 ||IS1[2]>LTRfragLen+4){
						#print "  " IS1[2], LTRfragLen-150 ,LTRfragLen,LTRfragLen+4
					#	print $0 "-> Deletion smallest than admitted (" LTRfragLen-150")"  > fileName[1]"_"ml".is.removed"
					#	next;
					#}
					
					
					if($2==16) {
						delta=IS1[1]
						sequenza=substr($10,1,length($10)-IS1[2])
					}	
					else {
						delta=0;
						sequenza=substr($10,IS1[2])
					}
					if(IS1[1]>=ml)
						print $1" "$3" "$4+delta" "$2" "$6" "IS1[1]" "sequenza> fileName[1]"_"ml".is.r1"
						#print $1" "$3" "$4+delta" "$2" "$6" "IS1[1]" "sequenza
					else {
						print $0 "-> Seq length smaller than " ml   > fileName[1]"_"ml".is.removed.r"
						next
					}
					
				}
				else if (length(cigar2)==2){
					#XSYM
					# If the 22 (+-4) bases are not found the IS is removed
		        		#if (cigar1[1]<LTRfragLen-150 ||cigar1[1]>LTRfragLen+4) {
					#	print $0 "-> Deletion smallest than admitted (" LTRfragLen-150")"  > fileName[1]"_"ml".is.removed"
					#	next;
					#}
					
					
					split(cigar1[2],IS2,"M")
					if($2==16) {
						delta=IS2[1]
						sequenza=substr($10,1,length($10)-cigar1[1])
					}	
					else {
						delta=0;
						sequenza=substr($10,cigar1[1])
					}
					if(IS2[1]>=ml)
						print $1" "$3" "$4+delta" "$2" "$6" "IS2[1]" "sequenza> fileName[1]"_"ml".is.r1"
						#print $1" "$3" "$4+delta" "$2" "$6" "IS2[1]" "sequenza
					else {
						print $0 "-> Seq length smaller than " ml   > fileName[1]"_"ml".is.removed.r"
					        next
					}
				}		

		}
		print $0 "-> No rules for this sequence"   > fileName[1]"_"ml".is.removed.r"
	}
}
