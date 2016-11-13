SCRIPTS=/home/extuser/bin/pIRS_111/
inputSet=$1
diploidSet=$2
lenFragment=$3
lenReads=$4
coverage=$5
stDev=$6
outputDir=$7
outName=${coverage}_${stDev}
echo "$SCRIPTS/pirs simulate -i $inputSet -I $diploidSet -m $lenFragment -l $lenReads -x $coverage -v $stDev -c 1 -o $outName > runInfo.o 2> runInfo.e"
$SCRIPTS/pirs simulate -i $inputSet -m $lenFragment -l $lenReads -x $coverage -v $stDev -c 1 -o $outName > runInfo.o 2> runInfo.e
