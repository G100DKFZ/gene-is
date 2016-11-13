outBase=$1
scriptDir=$2
lenFragment=300
primus=1

# set read length (for length > 100 you have to provide the corresponding error profile to pIRS)
for lenReads in 100 ;
do
	# set reads coverage
	#for coverage in 1 5 25 50 100 250 500 1000;
	for coverage in 250
	do
		cd $outBase/${lenFragment}_${lenReads}_${coverage}
		awk -v vector=AMT_Provirus -f $scriptDir/toBEDformat.awk resultsNoDup.csv |awk -v vector=AMT_Provirus -v header=$primus -v tabular=1 -v fileIn=$outBase/IntegrationSitesReal.txt -f $scriptDir/ISatitsics.awk >> $outBase/resultsStatistics_${lenFragment}_${lenReads}.csv
		for i in {1..150}; 
		do 
			awk -v thr=$i '{if ($4>=thr && $1!="AMT_Provirus") for (i=0;i<$4;i++) print $0}' resultsNoDup.csv.total.results|awk -v vector=AMT_Provirus -v header=0 -v tabular=1 -v fileIn=../IntegrationSitesReal.txt -f $scriptDir/ISatitsics.awk 
		done
		primus=0
	done
done
