# extractSCLIP.awk 2016
function resetForward() {
	nameF=""
	chrF=""
	startF=0
	leftF=0
	rigthF=0
	seqF=""
	qualF=""
	cigarF="100M"
}

function resetReverse() {
	nameR=""
	chrR=""
	startR=0
	leftR=0
	rigthR=0
	seqR=""
	qualR=""
	cigarR="100M"
}



function getSclip(flag,cigar) {
	split(cigar,vect,"S")
	# if the read contains more than 2 SOFTCLIPs discard it
	if(length(vect) > 3) {
		sclip[1]=0
		sclip[2]=0
		return 0
	}
	# the read contains 0, 1 or 2 softclip
	else {
		# the read's cigar has a perfect match: spilt(XM,S)
		if(length(vect)==1) {
			sclip[1]=0
			sclip[2]=0
			return 1
		} else {
			# the cigar has a left Sclip: split(YSXM,S) -> [Y] [XM]
			if(length(vect)==2 && vect[2]!="") {
				# Check the minimal length of the Sclip
				if(vect[1]>=thr)
					sclip[1]=vect[1]
				else
					sclip[1]=0
          sclip[2]=0

			# the read has a rigth Sclip: split(XMYS,S) -> [XMY] [] -> [X][Y]
			} else if(length(vect)==2 && vect[2]==""){

				split(vect[1],appoggio,"M")
				sclip[1]=0
				# extract the correct sclip dimesion
				if(appoggio[length(appoggio)]>=thr)
                                	sclip[2]=appoggio[length(appoggio)]
				else
					sclip[2]=0

      # the read has a left/rigth Sclip: split(ZSXMAIBMYS,S) -> [Z][XMAIBMY] [] -> [X][AIB][Y]
			} else {
				split(vect[2],appoggio,"M")
				if(vect[1]>=thr)
                                	sclip[1]=vect[1]
				else
					sclip[1]=0
				# extract the rigth sclip dimesion
				if (appoggio[length(appoggio)]>=thr)
          sclip[2]=appoggio[length(appoggio)]
				else
					sclip[2]=0
			}
			return 1
		}
	}
}


# return the not of a.
function not(       a) {
	a>0?a=1:a=0
	return 1-a
}

# This function implement the logic of well formatted Sclip.
# In our contest is a non matching region that belong to a pair end alignment
# that is in the correct position.

function printFragment(leftF,leftR,rigthF,rigthR,strand) {
	#print "printFragment -> logic"

  (leftF) ? a=1:a=0
	(leftR) ? c=1:c=0
	(rigthF) ? b=1:b=0
	(rigthR) ? d=1:d=0
###### gawk 4.X compatibility on ternary logical operators
	#strand=="+" ? valid=or(and(not(a),not(c),d),and(a,not(b),not(d))):valid=or(and(not(b),c,not(d)),and(not(a),b,not(c)))
###### gawk 3.x compatibility 
	strand=="+" ? valid=or(and(and(not(a),not(c)),d),and(and(a,not(b)),not(d))):valid=or(and(and(not(b),c),not(d)),and(and(not(a),b),not(c)))

	strand=="+" ? num_strand=1:num_strand=0
	if(valid==1) {
	    potentialReads++
	    #print "printFragment -> is valid",leftF,rigthF,leftR,rigthR
	    if (and(a,c)) {
		#print "printFragment -> both +"
        	if (leftF > leftR) {
                	chrom=chrF;start=startF;readName=nameF;sequence=seqF
			soft=substr(seqF,1,leftF);qual=substr(qualF,1,leftF);fragment=1
                } else {
                	chrom=chrR;start=startR;readName=nameR;sequence=seqR
			soft=substr(seqR,1,leftR);qual=substr(qualR,1,leftR);fragment=1
		}
            } else if (and(b,d)) {
		#print "printFragment -> both -"
                if (rigthF > rigthR) {
			seqLenF=length(seqF)
                	chrom=chrF;start=startF;readName=nameF;sequence=seqF
			soft=substr(seqF,seqLenF-rigthF);qual=substr(qualF,seqLenF-rigthF);fragment=2;
                        matchLenF=seqLenF-leftF-rigthF
                        start=startF+matchLenF
                } else {
			seqLenR=length(seqR)
                	chrom=chrR;start=startR;readName=nameR;sequence=seqR
			soft=substr(seqR,seqLenR-rigthR);qual=substr(qualR,seqLenR-rigthR);fragment=2;
                        matchLenR=seqLenR-leftR-rigthR
                        start=startR+matchLenR
		}
            } else if(a) {
			#print "printFragment -> a"
                	chrom=chrF;start=startF;readName=nameF;sequence=seqF
			soft=substr(seqF,1,leftF);qual=substr(qualF,1,leftF);fragment=1;
            } else if(d) {
			#print "printFragment -> d"
			seqLenR=length(seqR)
                	chrom=chrR;start=startR;readName=nameR;sequence=seqR
			soft=substr(seqR,seqLenR-rigthR);qual=substr(qualR,seqLenR-rigthR);fragment=2;
                        matchLenR=seqLenR-leftR-rigthR
                        start=startR+matchLenR
            } else if(c) {
			#print "printFragment -> c"
                	chrom=chrR;start=startR;readName=nameR;sequence=seqR
			soft=substr(seqR,1,leftR);qual=substr(qualR,1,leftR);fragment=1;
                        #start,stop,readLen,cigar,soft,qual,sequence=get_read_features(read2,firstM,0)
            } else if(b) {
			#print "printFragment -> b"
			seqLenF=length(seqF)
                	chrom=chrF;start=startF;readName=nameF;sequence=seqF
			soft=substr(seqF,seqLenF-rigthF);qual=substr(qualF,seqLenF-rigthF);fragment=2;
                        matchLenF=seqLenF-leftF-rigthF
                        start=startF+matchLenF
	    }
	    #print chrom "\t" start "\t" strand "\t" readName "\t" soft "\t" qual "\t" sequence "\t" num_strand "\t" fragment > prefix".sclip.txt"
	    print chrom "\t" start "\t" strand "\t" readName "\t" soft "\t" qual "\t" sequence "\t" num_strand "\t" fragment
	    return
	}
	return
}

BEGIN{
	# Setting the principal SAMTOOLS flags
	ALIGNED=1
	PROPER_PAIR=2
	THIS_REVERSE=16
	MATE_REVERSE=32
	NOT_PRIMARY=256
	FORWARD=64
	REVERSE=128

	PROPER_PAIR_NOT_PRIMARY=258
	totalReads=0
	properReads=0

  if(prefix=="")
    prefix="reads"
}
{
    	
    totalReads++
    if (specific==2) { 
	sclip[1]=0
	flag=$2; cigar=$6
	if(and(flag,PROPER_PAIR_NOT_PRIMARY)==2) {
		if (and(flag,FORWARD)) {
			if(and(flag,THIS_REVERSE))
				strand="-"
			else
				strand="+"
			success=getSclip(flag,cigar)
			nameF=$1
			chrF=$3
			startF=$4-1
			leftF=sclip[1]
			rigthF=sclip[2]
			seqF=$10
			qualF=$11
			cigarF=cigar
			reads+=1
			#forwSclip=sclip

		} else {
			if(and(flag,MATE_REVERSE))
				strand="-"
			else
				strand="+"

			success=getSclip(flag,cigar)
			nameR=$1
			chrR=$3
			startR=$4-1
			leftR=sclip[1]
			rigthR=sclip[2]
			seqR=$10
			qualR=$11
			cigarR=cigar
			reads+=1
			#revSclip=sclip
			#print $1,flag,cigar,sclip
		}
		if (reads==2 && nameF==nameR) {
			properReads++
			reads=0
			printFragment(leftF,leftR,rigthF,rigthR,strand)
		}

	}
    } else {
	sclip[1]=0
	flag=$2; cigar=$6
	if(and(flag,ALIGNED)) {
		if (and(flag,FORWARD)) {
			resetReverse()
			if(and(flag,THIS_REVERSE))
				strand="-"
			else
				strand="+"
			success=getSclip(flag,cigar)
			nameF=$1
			chrF=$3
			startF=$4-1
			leftF=sclip[1]
			rigthF=sclip[2]
			seqF=$10
			qualF=$11
			cigarF=cigar

		} else {
			resetForward()
			if(and(flag,MATE_REVERSE))
				strand="-"
			else
				strand="+"

			success=getSclip(flag,cigar)
			nameR=$1
			chrR=$3
			startR=$4-1
			leftR=sclip[1]
			rigthR=sclip[2]
			seqR=$10
			qualR=$11
			cigarR=cigar
		}
		# medium specificity. No regards on Alignment topology but reads should be in a proper pair
		if (specific==1 && and(flag,PROPER_PAIR)==0)
			next

		# low specificity. 
		properReads++
		printFragment(leftF,leftR,rigthF,rigthR,strand)
	}	
    }
}
END {
	print "Total aligned Reads: ",int(totalReads/2) 2>"stat.txt"
	print "Proper aligned Reads: ",properReads 2>"stat.txt"
	print "Potential vector-genome Reads: ",potentialReads 2>"stat.txt"
}
