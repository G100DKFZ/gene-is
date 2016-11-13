function evaluate(best,strand,posizione,prec,len,seq,name) {
	split(best,campi,"\t")
	if(strand=="+") {
		if(campi[9]>campi[10]) {
			genoPos=campi[9];
			genoStr="+"
		}
		else {
			genoPos=campi[10]
			genoStr"-"
		}
	} else {
		if(campi[9]>campi[10]) {
			genoPos=campi[10];
			genoStr="+"
		}
		else {
			genoPos=campi[9]
			genoStr="-"
		}
	}
	print prec,campi[2],name,genoPos,posizione,genoStr,strand,campi[3],len,seq
}
BEGIN {
	best=""
}
{
	split($1,vect,"@");
	if(prec!=vect[1]) {
		if(best!="") {
			evaluate(best,strand,posizione,prec,len,seq,name);
		}
		best=$0;
		score=$NR;
		prec=vect[1]
		name=vect[2]
		strand=vect[4]
		posizione=vect[3]
		len=vect[5]
		seq=vect[6]
	}
	else {
		if(best!="") {
			if ($NR==score)
				if ($2!=vectorStr)
					best="";
		}
	}
}
