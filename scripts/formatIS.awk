{
	if($2==vectorName && $3==vectorName) next
	if ($3==vectorName) chr[$2"\t"$4"\t"$7]++ 
	else chr[$3"\t"$5"\t"$6]++
}END{
	for (key in chr) print key"\t"chr[key]
}
