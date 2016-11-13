function evaluate(best,strand,posizione,id,len,seq,name,side,fragment) {
	split(best,campi,"\t")
	# the pairs were on the vector. So we have to detail IS position and orientation because the IS info is in this alignment. The orientation works only for LTR. In ITR is not reliable.
	if (name==vectorStr){
		# when fragment=2 the sclip is the rigth one
		# so the IS is the proximal. 
		if (fragment==2) {
			genoPos=campi[9]
			if(campi[9]>campi[10]) {
				if (1+side==1)
                        		genoStr="-"
				else
					genoStr="+"
			}
			else {
				if (0+side==1)
                        		genoStr="-"
				else
					genoStr="+"
			}
		# fragment=1 so the sclip is the left one
		# The IS is the disttal
                } else {
                        genoPos=campi[10]
			if(campi[9]>campi[10]) {
				if (1+side==1)
                        		genoStr="-"
				else
					genoStr="+"
                        }
                        else {
				if (0+side==1)
                        		genoStr="-"
				else
					genoStr="+"
                        }
		}
	}
	# the pairs are on the genome. The IS information is in the original alignment 
	else {
		if(side==0) {
			if(campi[9]>campi[10]) {
				genoPos=campi[10];
				genoStr="-"
			}
			else {
				genoPos=campi[10]
				genoStr="+"
			}
		} else {
			if(campi[9]>campi[10]) {
				genoPos=campi[9];
				genoStr="+"
			}
			else {
				genoPos=campi[9]
				genoStr="-"
			}
		}
	}
	print id,campi[2],name,genoPos,posizione,genoStr,strand,ostia,campi[3],len,seq
}
BEGIN {
	best=""

}
{
	split($1,vect,"@");

	if(prec!=$1) {
		if(best!="") {
			evaluate(best,strand,posizione,id,len,seq,name,side,fragment);
		}
		id=vect[1]
		best=$0;
		score=$NR;
		prec=$1
		name=vect[2]
		strand=vect[4]
		posizione=vect[3]
		len=vect[5]
		seq=vect[6]
		side=vect[7]
		fragment=vect[8]
	}
	else {
		if(best!="") {
			if ($NR==score)
				# most of the vector aligned sequeces are found in both ITR. So same score in vector is allowed
				if ($2!=vectorStr)
					best="";
		}
	}
}
END  {
	if(best!="") {
		evaluate(best,strand,posizione,id,len,seq,name,side);
	}
}	
