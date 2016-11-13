function average(vector) {
	count=0
	mu=0
	for (keys in vector) {
		mu+=vector[keys]
		#count++
	}
	#mu=mu/count

	return mu/(length(vector)+0.0000000001)
}

function stdev(vector) {
	count=0
	mu=0
	st=0
	mu=average(vector)
	for (keys in vector) {
		if(vector[keys]>0)
			st+=(vector[keys]-mu)*(vector[keys]-mu)
		#count++
	}
	#st=sqrt(st/(count-1.0000000))
	return sqrt(st/(length(vector)-0.9999999999))
}

function red(s) {
	    return "\033[1;31m" s "\033[0m "
}

#function white(s) {



BEGIN{
	# strandError contains the number of IS that are annotated with conflicting strand
	strandError=0
	while(getline < fileIn >0) {
		IS[$2]=0
		isID[$2]=$1" "$2
	}
}
{
	for (keys in IS) {
		#print "searching position: " keys " - "$2 " = " $2-keys
		if(($2-keys)*($2-keys)<10000) {
			# if verbose==2 show also the found IS True Positives (TP)
			if(verbose>1) print $1,$2,"TP"
			stat[keys" "$2]++
			mean[keys]+=$2-keys
			IS[keys]++
			totalTrueReads++
			found=1
			break
		}
	}
	if (!found) {
	        if(verbose>=1) print $1,$2,"FP"	
		False[$2]++;
		totalFalseReads++
	}
	if(match($0,"strand_error")) strandError++
	found=0
		
}
END{
	miss=0
	found=0
	for (keys in IS) {
		if (IS[keys]>0) found++;
		else if (keys >1) {
			miss++;
			if(verbose>=1) print isID[keys],"FN"
		}
	}
	total=found+miss
	for (keys in False)
		falseCount++


	for (keys in mean) {
		mean[keys]=mean[keys]/IS[keys]
		#print "count " keys " " IS[keys]
	}
	for (realIS in mean) {
		for (keys in stat) {
			split(keys,vect," ")
			if (vect[1]==realIS){
				if(max[realIS]<sqrt((vect[2]-realIS)*(vect[2]-realIS)))
				   max[realIS]=sqrt((vect[2]-realIS)*(vect[2]-realIS))
				#std[realIS]+=((vect[2]-realIS)-mean[realIS])*((vect[2]-realIS)-mean[realIS])
				#print realIS " - " vect[1] " - " vect[2] " - " realIS " "mean[realIS] " - " count[realIS]
			}
		}
		#print realIS" -> " count[realIS], mean[realIS], std[realIS]
		#print "IS: "realIS "\treads: "IS[realIS]"\tmean: "mean[realIS]"\tvariance: " sqrt(std[realIS]/(IS[realIS]-1.00001))
		if (massimo<max[realIS])
			massimo=max[realIS]
	}	
	if (!tabular) {
		print "*** Summary ***"
		print "Total number of reads: "red(totalTrueReads+totalFalseReads)
		print "Number of reads in true IS: "red(totalTrueReads)
		print "Average number of reads in true IS (stdev): " red(average(IS)) " ("red(stdev(IS))")"
		print "Number of reads in false IS: "red(totalFalseReads)
		print "Average number of reads in false IS (stdev): " red(average(False)) " ("red(stdev(False))")"
		print "Reads Signal/Noise:" red(totalTrueReads/(totalFalseReads+0.01))
		print "Sensibility :" red(totalTrueReads/(totalTrueReads+totalFalseReads))
		print "Total IS: " red(found+miss)
		print "Found " red(found) " true IS (" red(found/total*100) "; IS Sensibility)"
		print "Missed " red(miss) " true IS (" red(miss/total*100) ")"
		print "Number of False positive IS: " red(falseCount)
		print "Average maximal distance from Precise IS (std): "red(average(max))" ("red(stdev(max))")"
		print "Maximal found distance: " red(massimo)
		print "Number of conflicting IS strands: " red(strandError)
	}	
	else {
		if (header) {
			print "#reads_with_IS\t#reads_in_true_IS\taverage\tstd\t#reads_in_false_IS\taverage\tstd\tS/N\tSensibility\t#IS\tFound\tPercentage\tMissed\tPercentage\t#False_Positive\taverage_cluster_distance\tstd\tmaximal\tstrandError"
			print totalTrueReads+totalFalseReads"\t"totalTrueReads"\t" average(IS)"\t"stdev(IS)"\t"totalFalseReads"\t"average(False)"\t"stdev(False)"\t"totalTrueReads/(totalFalseReads+0.01)"\t"totalTrueReads/(totalTrueReads+totalFalseReads)"\t"found+miss"\t"found"\t"found/total*100"\t"miss"\t"miss/total*100"\t"falseCount"\t"average(max)"\t"stdev(max)"\t"massimo"\t"strandError
		}
		else
			print totalTrueReads+totalFalseReads"\t"totalTrueReads"\t" average(IS)"\t"stdev(IS)"\t"totalFalseReads"\t"average(False)"\t"stdev(False)"\t"totalTrueReads/(totalFalseReads+0.01)"\t"totalTrueReads/(totalTrueReads+totalFalseReads)"\t"found+miss"\t"found"\t"found/total*100"\t"miss"\t"miss/total*100"\t"falseCount"\t"average(max)"\t"stdev(max)"\t"massimo"\t"strandError

	}
}
