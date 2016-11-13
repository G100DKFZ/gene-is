# SortLinkerParallel 1.0
# 
# Sorting of Miseq run giving a bouple barcodie key found in tagtable file. The input raw seq files should be asymetric (450/50)
# Output could be fasta or fastq (variable fa)

BEGIN{ 
	format = "%a %b %e %H:%M:%S %Z %Y"
	bar=0
	stretchA=""
	stretchQ=""
	
	if (outDir=="")
		outDir="."

	print strftime(format) >> outDir"/sort.log"
	while ((getline line <fileIn)>0){
		tagTab[line]=1
	}
	close(fileIn)
	print fileIn
		
#	FS="\n";RS="\n@HWI";
	FS="\n";RS=readHeader;
	print RS
}


FNR==NR {
	split($1,ident," ");    # the first filed (line) has the reads name. Used to index the hash table
	seq[ident[1]]=$2;	# the sequence is on line 2
	qual[ident[1]]=$4	# the quality on line 4
	next
}

{

	#print "Uguali" 
	split($1,ident," ") 	# index the hash table with the read id of the other read
	flag=0
	for (line in tagTab) {	# for each couple of barcodes
                split(line, vect,"\t") 	# split the barcodes stucture
		#print substr(seq[ident[1]],1,20),vect[1]"\t"$2,vect[2]
                if(vect[2]==""){	# if no barcode on field 2 of tag table
			 if (match(substr(seq[ident[1]],1,20),vect[1])) {
	                         goodR1=seq[ident[1]]
				 if(!fa) {
	                         	goodR1qual=substr(qual[ident[1]],1,match(seq[ident[1]],vect[3])-1)
	                         	print "@HWI" ident[1]"\n"goodR1"\n+\n"goodR1qual >> outDir "/"vect[4]"__" procID #| gzip
				 }
				 else 
					 print ">HWI" ident[1]"\n"goodR1 >> outDir "/"vect[4]"__" procID #| gzip
	                         counterFiles[vect[4]]++;
	                         print "ok -> "$1 >> outDir "/sort__" procID
	                         flag=1
	                         break
	                 }
	        }
		else if ((match(substr(seq[ident[1]],1,20),vect[1]))&&(match($2,vect[2]))){      	# if in the first 20 characters of the forward the LTR barcode is present and also the lnk barcode 
													# is present on the reverse 
				if (trimming) { 	# if the flag trimming is true then this script trimm the sequence from the barcote to the end of the sequence
	                        	goodR1=substr(seq[ident[1]],1,match(seq[ident[1]],vect[3])-1)		# remove the DNA string that exceed the linker cassette
	                        	goodR1qual=substr(qual[ident[1]],1,match(seq[ident[1]],vect[3])-1)	# also on the quality string
				} else {
					goodR1=seq[ident[1]]
					goodR1qual=qual[ident[1]]
				}
				lenR1=length(goodR1)
	                        if ((lenR1<=100)&&(lenR1>=55)){
	                                goodR1=goodR1 stretchA
					if(!fa) 
	                                	goodR1qual=goodR1qual stretchQ
	                        }
		                if (lenR1<55) {
		                        goodR1=seq[ident[1]]
					if(!fa) 
		                        	goodR1qual=qual[ident[1]]
	                        }
				if(!fa) 
	                        	print "@HWI" ident[1]"\n"goodR1"\n+\n"goodR1qual >> outDir "/"vect[4]"__" procID #| gzip
				else 
	                        	print ">HWI" ident[1]"\n"goodR1 >> outDir "/"vect[4]"__" procID #| gzip
		                counterFiles[vect[4]]++;
		                print "ok -> "$1 >> outDir "/sort__" procID
		                flag=1
		                break
		}

    	}
	if (flag==0){
                leftL=substr(seq[ident[1]],1,15)
                lnkPos=match($2,"CCTAAC")
	    	leftR=substr($2,1,15)	
                if(lnkPos>=10)
                        leftR=substr($2,lnkPos-10,15)
                leftOver[leftL"\t"leftR]++
                print "Left-> "leftL,leftR,$1 >> outDir"/sort__" procID
                print "@HWI" ident[1]"\t"seq[ident[1]] >>outDir "/leftover__" procID
                count ++
        }
}
END{
	print strftime(format) >> outDir "/sort__" procID
	for (wrongBC in leftOver)
	    if(leftOver[wrongBC]>10)
	            print wrongBC "\t" leftOver[wrongBC] >> outDir "/wrongBC__" procID
}
