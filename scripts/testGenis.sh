## GENE-IS testing script
## Sept 2016

function return_error
{
	echo "!!! Assertion failed !!!"
	echo "$@" 1>&2
	echo $COMMENT_LINE
	exit 1
}

OUT_TES_PAIR_FILE="$GENIS/test/targetedSequencing/results/pairedEnd/testDataTS.ResultsClusteredAnnotated.csv"
OUT_TES_SINGLE_FILE="$GENIS/test/targetedSequencing/results/singleEnd/testDataTS.ResultsClusteredAnnotated.csv"
TEST_TES_PAIR_FILE="$GENIS/test/targetedSequencing/results/testDataTS.pair.csv"
TEST_TES_SINGLE_FILE="$GENIS/test/targetedSequencing/results/testDataTS.single.csv"
OUT_LAM_PAIR_FILE="$GENIS/test/LAM-PCR/results/testRunLAM.25-0.95-0.9-3-1000.ResultsClusteredAnnotated.csv"
TEST_LAM_PAIR_FILE="$GENIS/test/LAM-PCR/testDataLAM.pair.csv"


ERROR_ASSERT_NO_GENIS="Set the environment variable $GENIS (i.e. export GENIS=/path_to_/GENE-IS)"
ERROR_ASSERT_TES_PAIR_DIFF="Output File ${TEST_TES_PAIR_FILE} is not the same as expected"
ERROR_ASSERT_TES_SINGLE_DIFF="Output File ${TEST_TES_SINGLE_FILE} is not the same as expected"
ERROR_ASSERT_LAM_PAIR_DIFF="Output File ${TEST_LAM_PAIR_FILE} is not the same as expected"
ERROR_ASSERT_TES_NO_FILE="Output file ${OUT_TES_PAIR_FILE} not existing"
ERROR_ASSERT_LAM_NO_FILE="Output file ${OUT_LAM_PAIR_FILE} not existing"
COMMENT_LINE="###################################################################"

if [ -z "$GENIS" ]; then
	echo $COMMENT_LINE
	return_error $ERROR_ASSERT_NO_GENIS
	echo $COMMENT_LINE
fi

PS3='Please enter your choice: '
options=("Targeted Sequencing Pair BWA" "Targeted Sequencing Single" "LAM-PCR" "All" "Clear" "Quit")
select opt in "${options[@]}"
do
	case $opt in
  #         "Targeted Sequencing Pair Bowtie2")
#		rm $GENIS/test/targetedSequencing/results/pairedEnd/*
 #               perl -I $GENIS/lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/targetedSequencing/results/pairedEnd/ -c $GENIS/configFileTestAgilentPairBowtie2.txt 
                #perl -I $GENIS/lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/AGILENT/results/pair/ -c $GENIS/configFileTestAgilentPair.txt &>> $GENIS/test/AGILENT/results/pair/log.txt 
#            	;;
           "Targeted Sequencing Pair BWA")
		rm $GENIS/test/targetedSequencing/results/pairedEnd/* 2>/dev/null
                perl -I $GENIS/lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/targetedSequencing/results/pairedEnd -c $GENIS/configFile_targetedSequencing_pairedEnd.txt
		echo "##################### Testing Pair-ends Output #########################"
		echo "Generated Output File $TEST_TES_PAIR_FILE"
		echo "Template Output File $OUT_TES_PAIR_FILE"
		if [ ! -f $OUT_TES_PAIR_FILE ] ; then
			return_error "$ERROR_ASSERT_TES_NO_FILE"
		fi
		if  ! cmp $TEST_TES_PAIR_FILE   $OUT_TES_PAIR_FILE > /dev/null 2>/dev/null; then 
			return_error "$ERROR_ASSERT_TES_PAIR_DIFF"
		fi
		echo $COMMENT_LINE
		echo "Targeted Sequencing Pair worked as expected!"
		echo $COMMENT_LINE
                #perl -I $GENIS/lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/AGILENT/results/pair/ -c $GENIS/configFileTestAgilentPair.txt &>> $GENIS/test/AGILENT/results/pair/log.txt 
            	;;
           "Targeted Sequencing Single")
		rm $GENIS/test/targetedSequencing/results/singleEnd/* 2>/dev/null
                perl -I $GENIS/lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/targetedSequencing/results/singleEnd/ -c $GENIS/configFile_targetedSequencing_singleEnd.txt
		echo "##################### Testing Single-end Output #########################"
		echo "Generated Output File $TEST_TES_SINGLE_FILE"
		echo "Template Output File $OUT_TES_SINGLE_FILE"
		if [ ! -f $OUT_TES_single_FILE ] ; then
			return_error "$ERROR_ASSERT_TES_NO_FILE"
		fi
		if  ! cmp $TEST_TES_SINGLE_FILE   $OUT_TES_SINGLE_FILE > /dev/null 2>/dev/null; then 
			return_error "$ERROR_ASSERT_TES_SINGLE_DIFF"
		fi
		echo $COMMENT_LINE
		echo "Targeted Sequencing Single-end  worked as expected!"
		echo $COMMENT_LINE
            	;;
           "LAM-PCR")
		rm -r $GENIS/test/LAM-PCR/results/* 2>/dev/null
               	perl -I $GENIS/lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/LAM-PCR/results -c $GENIS/configFile_LAM-PCR_pairedEnd.txt
                echo "##################### Testing LAM-PCR Pair-ends Output #########################"
                echo "Generated Output File $TEST_LAM_PAIR_FILE"
                echo "Template Output File $OUT_LAM_PAIR_FILE"
                if [ ! -f $OUT_LAM_PAIR_FILE ] ; then
                        return_error "$ERROR_ASSERT_LAM_NO_FILE"
                fi
                if  ! cmp $TEST_LAM_PAIR_FILE   $OUT_LAM_PAIR_FILE > /dev/null 2>/dev/null; then
                        return_error "$ERROR_ASSERT_LAM_PAIR_DIFF"
                fi
                echo $COMMENT_LINE
                echo "LAM-PCR Pair worked as expected!"
                echo $COMMENT_LINE
		;;


	   "All")
		rm $GENIS/test/targetedSequencing/results/pairedEnd/*
                perl -I ../lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/targetedSequencing/results/pairedEnd/ -c $GENIS/configFile_targetedSequencing_pairedEnd.txt
		rm $GENIS/test/targetedSequencing/results/singleEnd/*
                perl -I ../lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/targetedSequencing/results/singleEnd/ -c $GENIS/configFile_targetedSequencing_singleEnd.txt
		rm $GENIS/test/LAM/results/*
               	perl -I $GENIS/lib/ $GENIS/scripts/GENIS.pl -o $GENIS/test/LAM-PCR/results -c $GENIS/configFile_LAM-PCR_pairedEnd.txt
		;;
	   "Clear")
		rm $GENIS/test/targetedSequencing/results/pairedEnd/* 2>/dev/null
		rm $GENIS/test/targetedSequencing/results/singleEnd/* 2>/dev/null
		rm -r $GENIS/test/LAM-PCR/results/* 2>/dev/null
		;;
           "Quit")
                break
            	;;
           *) 
		echo Wrong option [1-7]
		;;
        esac
done
