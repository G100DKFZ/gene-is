# synthExperiment.sh
#
# raffaele.fronza@nct-heidelberg.de
# 09.2014
# Description:
# The script simulate a complete ngs analysis using GENIS on a synthetic set of DNA fragments (that contains known IS)
# Input args:
# outputBase = output base directory
# inputSet = DNA fragments used to simulate illumina reads
# diploidSet = Diploid version of inputSet (i.e. generated via pIRS)
# configTemplateFile = template file for GENIS pipeline
# scriptDir = base directory of GENIS scripts

# Dependencies:
# A functional version of GENIS.
# The script synthReds.sh for the generation of synthetic illumina reads
# toBEDformat.awk that translate resultsNoDup.csv (raw found IS set) to BED format (chr, pos, pos+1, segqence, %identity, strand)
# IStatistics.awk perform some statistics on tabular files in (2-col) BAD file (che,pos)

# Usage example
# bash ../scripts/synthExperiment.sh ~/genis/testSynth ~/genis/testSynth/RandomData_100frags_10000bp_appendedWithoutN.fasta ~/genis/testSynth/HumanpVD43.snp.indel.inversion.fa   ~/genis/testSynth/configFileTestAgilentPairBWASynth.txt /home/extuser/genis/scripts  ~/genis/testSynth/IntegrationRealRepetitive.txt
#

outBase=$1
inputSet=$2
diploidSet=$3
configTemplateFile=$4
scriptDir=$5
realSet=$6
lenFragment=200
primus=1

# set read length (for length > 100 you have to provide the corresponding error profile to pIRS)
for lenReads in 100 ;
do
	# set reads coverage
	for coverage in  1 10 100 1000 ;
 	do
		# check existence of output directory. Remove content if present
		if [ -d "$outBase/${lenFragment}_${lenReads}_${coverage}" ]; then
			echo "The output dir exists. Removing the content."
			echo "rm $outBase/${lenFragment}_${lenReads}_${coverage}/*"
			rm  $outBase/${lenFragment}_${lenReads}_${coverage}/*
			  #exit
		fi
		# Perform the synthetic Illumina run with the imposed fregment, read length and coverage 
		echo "mkdir $outBase/${lenFragment}_${lenReads}_${coverage}"
		mkdir $outBase/${lenFragment}_${lenReads}_${coverage}
		echo "cd $outBase/${lenFragment}_${lenReads}_${coverage}"
		cd $outBase/${lenFragment}_${lenReads}_${coverage}
		echo "bash $scriptDir/synthReads.sh $inputSet $diploidSet $lenFragment $lenReads $coverage 10 ${coverage}_10"
		bash $scriptDir/synthReads.sh $inputSet $diploidSet $lenFragment $lenReads $coverage 10 ${coverage}_10

		# Modify the selected GENIS configuration file
		#$sedDir="$outBase/${coverage}_10_${lenReads}_${lenFragment}"
		echo "sed  -i \"s/\(forward *= *\).*/\1${coverage}_10_${lenReads}_${lenFragment}_1.fq.gz/\" $configTemplateFile"
		sed -i "s/\(forward *= *\).*/\1${coverage}_10_${lenReads}_${lenFragment}_1.fq.gz/" $configTemplateFile
		echo "sed  -i \"s/\(reverse *= *\).*/\1${coverage}_10_${lenReads}_${lenFragment}_2.fq.gz/\" $configTemplateFile"
		sed -i "s/\(reverse *= *\).*/\1${coverage}_10_${lenReads}_${lenFragment}_2.fq.gz/" $configTemplateFile

		# Run GENIS and detect IS
		echo "perl -I /home/saira2/GENE-IS_edited/lib $scriptDir/GENIS.pl -c $configTemplateFile -o $outBase/${lenFragment}_${lenReads}_${coverage}"
		perl -I /home/saira2/GENE-IS/lib $scriptDir/GENIS.pl -c $configTemplateFile -o $outBase/${lenFragment}_${lenReads}_${coverage}
		
		# Analyze results and write them in tabular format on $outBase/resultsStatistics_${lenFragment}_${lenReads}.csv 
		echo "awk -v vector=AMT_Provirus -f $scriptDir/toBEDformat.awk resultsNoDup.csv |awk -v vector=AMT_Provirus -v header=$primus -v tabular=1 -v fileIn=$realSet -f $scriptDir/ISatitsics.awk >> $outBase/resultsStatistics_${lenFragment}_${lenReads}.csv"
		awk -v vector=PEX-A2 -f $scriptDir/toBEDformat.awk resultsNoDupSingle.csv |awk -v vector=PEX-A2 -v header=$primus -v tabular=1 -v fileIn=$realSet -f $scriptDir/ISatitsics.awk >> $outBase/resultsStatistics_${lenFragment}_${lenReads}.csv
		primus=0
	done
done
