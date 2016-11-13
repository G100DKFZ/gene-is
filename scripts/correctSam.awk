# Prototype for correction of sam files

{
	#remove secondary alignment
	if($2 > 255 && NF>3)
		next
	# Uncomment to relax stringency on pairs distance
	#if(($2!=163 && $2!=99 && $2!=147 && $2!=83 && $2!=161 && $2!=97 && $2!=145 && $2!=81) && NF>3)
	
	# Uncomment to increase stringency on pairs distance
	if(($2!=163 && $2!=99 && $2!=147 && $2!=83) && NF>3)
		next
	print $0
}
