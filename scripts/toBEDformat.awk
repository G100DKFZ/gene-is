# This script is a prototype for rearrangemnts of results output to the BED format
# the basically rearrange the columns in order to have 
# 1 chromosome
# 2 position start
# 3 position end
# 4 reads sequence (full in this case)
# 5 identity score
# 6 strand

{
	if(vector==$2) 
		print "chr" $3,$5,$5+1,$NF,$8,$6
	else 
		print "chr" $2,$4,$4+1,$NF,$8,$7
}
