{
	if ((($2-IS)*($2-IS)<=18*18)&&(chrPrec==$1)) {
		if ($3==strandPrec) {
			if($3=="+") {
				if($4<count/10) 
					count+=$4
				else {
					print chrPrec,IS,"+",count,"*"
					IS=$2
					count=$4
				}
			}
			else if($3=="-") {
				if($4>precCount*10) {
					count+=$4
					IS=$2
					precCount=$4
				}
				else {
					count+=$4
				}
			}
		}
		else {
			print chrPrec,IS,strandPrec,count
			IS=$2
			count=$4
			strandPrec=$3
		}
	}
	else {
		if(first==1) 
		print chrPrec,IS,strandPrec,count
		first=1
		chrPrec=$1
	        IS=$2
		count=$4
		strandPrec=$3
	}
}
END {
	print chrPrec,IS,strandPrec,count
}
