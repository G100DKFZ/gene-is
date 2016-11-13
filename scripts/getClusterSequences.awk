{
	if(pos==$3) {
		if ($1>max) {
			max=$1; 
			out=$1
			seq=$6
		} else if ($1==max)
			out+=$1
	

	} else {
		print chr"_"pos,max,out,seq
		chr=$2;
		pos=$3;
		max=$1;
		out=$1
		seq=$6
	}
}
END{
	print chr"_"pos,max,out,seq
}
