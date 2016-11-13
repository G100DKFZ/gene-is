# filterSClipped.awk	2014
# Prototype 1.0
# Given a fastq file generated via extractSClip.pl it return all the reads that are longer than len (default 20)
# Files sorted by alignment position should be used
# in a file that is named in "fileName" variabile (default filtered)
# If "chr" is given the output file is a composition of "outfile chr"
# outDir is the working directory ="./" by default
BEGIN {
		
	if (!outDir)
		outDir="."
	if (!len)
		len=20
}

{
	if(length($5)>=len) { # length of the sequence match to the filter
		if ((chr)&&(chr==$1)) {	# only sequences found in chr are considered
			#print $0
			if (nodup) {			# remove identical sequences
				if (prec!=$5) {		# found adiacent in the list
					if (!chr){	# use default filename 
						if (fa)		# output in fasta format
							print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5 > outDir"/filtered" $1 ".fa"
						else 	# output in fastq format
							print "@"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5"\n+\n"$6 > outDir"/filtered" $1 ".fq"
					} else {
						#print $0	
						if (fa)
							print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5 > outDir"/"chr ".fa"
						else
							print "@"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5"\n+\n"$6 > outDir"/"chr ".fq"
					}
				}
			}
			else	{ # do adiacent duplications are considered
					if (!chr) {
						if (fa) 
							print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5 > outDir"/filtered" $1 ".fa"
						else
							print "@"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5"\n+\n"$6 > outDir"/filtered" $1 ".fq"
					} else	{
						if (fa)
							print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5 > outDir"/"chr $1 ".fa"
						else 
							print "@"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5"\n+\n"$6 > outDir"/"chr ".fq"
					}
			}
		} else if (!chr) { # all the is are considered
	       		if (nodup) { 	
				if (prec!=$5) {
					if(!fileName) {
						if (fa)		# output in fasta format
							print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5 > outDir"/filtered" ".fa"
						else 	# output in fastq format
							print "@"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5"\n+\n"$6 > outDir"/filtered" ".fq"
					} else	{
						if (fa)
							print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5 > outDir"/"fileName ".fa"
						else
							print "@"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5"\n+\n"$6 > outDir"/"fileName ".fq"
					}
				}
			}
			else	{# do adiacent duplications are considered
					if (!fileName) {
						if (fa) 
							print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5 > outDir"/filtered" ".fa"
						else
							print "@"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5"\n+\n"$6 > outDir"/filtered" ".fq"
					} else	{
						if (fa)
							print ">"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5 > outDir"/"fileName ".fa"
						else 
							print "@"$4"@"$1"@"$2"@"$3"@"length($5)"@"$7"@"$8"\n"$5"\n+\n"$6 > outDir"/"fileName ".fq"
					}
			}

		}
	}	
	prec=$5
}
