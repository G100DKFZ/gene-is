#!/usr/bin/perl -w
###################################################################################
#
#  GENE-IS, a fully automated pipeline for detection of viral integration sites
#  in host genomes in LAM-PCR or Targeted Sequencing based data, generated as a 
#  single or paired end reads
#           
#############################################################################
#
#  GENE-IS is free software
#
#  Version 1.0
#  Last update: 27/01/2016
#
#  GENIS.pl is part of GENE-IS. It is the interface of the software.
#
##############################################################################

use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use GeneisLib;

my @usage;
push @usage, "Program: GENE-IS, a tool for characterizing vector-genome insertion sites .\n";
push @usage, "Version: 1 \n\n";
push @usage, "Usage: GENIS.pl -c <configuration file> [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this information\n";
push @usage, "  -c, --config   Configuration file <required>\n";
push @usage, "  -o, --output   The directory to store software output, default is current working directory\n";


my $help;
my $config_file;
my $output_dir;
my $threads;
my $type;
my $genomeVector;
GetOptions
(
 'h|help|?'    => \$help,
 'config=s'    => \$config_file,
 'output=s'    => \$output_dir,
);


if ($help) {
   print @usage;
   exit(0);
}
if (defined $output_dir) {
    if (!-e $output_dir){
    	print "\nThe output directory $output_dir does not exist!\n\n";
        print @usage;
	exit;
    }
}else{
    $output_dir = getcwd;
}
if (defined $config_file) {
   if (!-e $config_file){
   	print "\nThe configuration file $config_file does not exist!\n\n";
        print @usage;
	exit;
   }
}else{
    print "Please provide a configuration file!\n\n";
    print @usage;
    exit;
}


#####################################################################################################

my $config = new();
$config->read($config_file);

my $paired=0;
my $lenLam="";
my $lamType="";
system "(echo  ......................................................)";
my $scriptDir = $config->get_value("scriptDir");
my $libDir=$config->get_value("libDir");

if ($config->has_value("type")) {
	$type=$config->get_value("type");
	if ( $type ne "AGILENT" ){
		$type="LAM";
	}
} else {
	$type="LAM";
}

system "date";
print "###################################################################\n";
print "Pre-proceesing (Quality Filtering and Adapter Trimming) in progress...\n\n";


my $forward = $config->get_value("forward");
my $samtools =  $config->get_value("samtools");
if ($type eq "LAM") {
	if ($config->has_value("lamType")){
		$lamType=$config->get_value("lamType");
		if ( $lamType eq "sorting") {
			$lamType="sorting";
		} 
		elsif ($lamType eq "extracting" )  {
			$lamType="extracting"
		}
		else {
			$lamType="both";			
		}
	}
	my $reverse = $config->get_value("reverse");
	my $tagTable = $config->get_value("tagTableName");
	my $group= $config->get_value("group");
	my $linker = $config->get_value("linker");
	my $megaprimer= $config->get_value("megaprimer");
	my $anchor = "noAnchor";
	my $anchorMM = 0.1;
	if ($config->has_value("anchor")){
		$anchor=$config->get_value("anchor");
		if ($config->has_value("anchorMM")){
			$anchorMM=$config->get_value("anchorMM");
		}
	}
	$lenLam= $config->get_value("len");
	my $barThr= $config->get_value("barThreads");
	my $lamPrefix= $config->get_value("lamPrefix");
	my $readHeader= $config->get_value("readHeader");
	my $qual = $config->get_value("qual");
	my $skewer =  $config->get_value("skewer");
	print "============================>>>> $lamType\n\n";
	system("echo \"perl $scriptDir/sortBarcode.pl -f $forward -r $reverse -tt $tagTable -o $output_dir -g $group -s $scriptDir -p $barThr -l $lenLam -sOut $lamPrefix -ty $lamType -anchor $anchor -anchorMM $anchorMM   -readHeader $readHeader -qual $qual -sk $skewer \"  ");
	#system("perl $scriptDir/sortBarcode.pl -f $forward -r $reverse -tt $tagTable -o $output_dir -g $group -s $scriptDir -p $barThr -lnk $linker -mega $megaprimer -l $lenLam -sOut $lamPrefix -ty $lamType");
	system("perl $scriptDir/sortBarcode.pl -f $forward -r $reverse -tt $tagTable -o $output_dir -g $group -s $scriptDir -p $barThr -l $lenLam -sOut $lamPrefix -ty $lamType -anchor $anchor -anchorMM $anchorMM   -readHeader $readHeader -qual $qual -sk $skewer " );

} else {
	my $qual = $config->get_value("qual");
#	my $perc = $config->get_value("perc"); 
	my $adaptF = $config->get_value("adaptF"); 
	my $adaptR = $config->get_value("adaptR"); 
	my $suffOut = $config->get_value("suffOut"); 
	my $skewer =  $config->get_value("skewer");
	my $reverse;
	if ($config->has_value("reverse")){
		$paired=1;
	        $reverse = $config->get_value("reverse");
	        system "(echo \"perl $scriptDir/filteringTrimming.pl  -f $forward -r $reverse -qual $qual  -adaptF $adaptF  -adaptR $adaptR -sOut $suffOut -o $output_dir -sk $skewer\"  )";
	        system "(perl $scriptDir/filteringTrimming.pl  -f $forward -r $reverse -qual $qual -adaptF $adaptF  -adaptR $adaptR -sOut $suffOut -o $output_dir -sk $skewer )";
	
	}else{
		system "(echo \"perl $scriptDir/filteringTrimming.pl -f $forward -qual $qual  -adaptF $adaptF  -sOut $suffOut  -o $output_dir  -sk $skewer\"  )";
		system "(perl $scriptDir/filteringTrimming.pl -f $forward -qual $qual  -adaptF $adaptF  -sOut $suffOut  -o $output_dir  -sk $skewer )";
	}
}

########################################################################################################
#print "###################################################################\n";
#print "Step 2 Indexing ...\n\n";

if ($config->has_value("index")){
        my $genome = $config->get_value("genome");
        my $vector= $config->get_value("vector");
        my $aligner = $config->get_value("aligner");
        my $indexOut = $config->get_value("indexOut"); 
	print "perl $scriptDir/index.pl -g $genome -v $vector  -a $aligner  -iOut $indexOut  -o $output_dir\n";
        system "(perl $scriptDir/index.pl -g $genome -v $vector  -a $aligner  -iOut $indexOut  -o $output_dir)";
	$genomeVectorIndex = "$output_dir/$indexOut"
}else{
#        print "You decided to use already available index files.\n";
	$genomeVectorIndex = $config->get_value("genomeVectorIndex");
}

########################################################################################################
########################################################################################################

print "\n###################################################################\n";
print "Alignment in process...\n\n";

if ($forwardOut) {
	my $forward=$forwardOut
}

if(!$config->has_value("threads")) { 
	$threads=1;
} else { 
	$threads = $config->get_value("threads");
}

# my $genomeVector = $config->get_value("genomeVectorIndex");
my $aligner = $config->get_value("aligner");
my $alignmentOut = $config->get_value("alignmentOut"); 
if ($type eq "AGILENT") {
	my $suffOut = $config->get_value("suffOut");
	print "=======>>>>>>>>>>>> $paired\n";
	if ($paired){
		if($reverseOut) {
			my $reverse=$reverseOut;
		}
       	 	system "(echo \"perl -I $libDir $scriptDir/alignment.pl  -p $threads -f $suffOut-pair1.fastq -r $suffOut-pair2.fastq  -gv $genomeVectorIndex  -a $aligner  -aOut $alignmentOut  -o $output_dir -t $type -s $scriptDir -sam $samtools \" )";
        	system "(perl -I $libDir $scriptDir/alignment.pl  -p $threads -f $suffOut-pair1.fastq -r $suffOut-pair2.fastq  -gv $genomeVectorIndex  -a $aligner  -aOut $alignmentOut  -o $output_dir -t $type -s $scriptDir -sam $samtools )";
	}else {
		system "(echo \"perl -I $libDir $scriptDir/alignment.pl  -p $threads -f $suffOut.fastq  -gv $genomeVectorIndex   -a $aligner  -aOut $alignmentOut  -o $output_dir -t $type -sam $samtools \" )";
		system "(perl -I $libDir $scriptDir/alignment.pl  -p $threads -f $suffOut.fastq  -gv $genomeVectorIndex   -a $aligner  -aOut $alignmentOut  -o $output_dir -t $type -sam $samtools )";
	}
} else {
	if ($lamType eq "both" or $lamType eq "extracting") {
		system "(echo \"perl -I $libDir $scriptDir/alignment.pl  -p $threads  -gv $genomeVectorIndex   -a $aligner -o $output_dir -s $scriptDir -t $type -sam $samtools \" )";
		system "(perl -I $libDir $scriptDir/alignment.pl  -p $threads  -gv $genomeVectorIndex   -a $aligner  -o $output_dir -s $scriptDir -t $type -sam $samtools )";
	}
	else {
		print "No Alignment step required\n"
	}
}	

########################################################################################################
########################################################################################################

print "###################################################################\n";
print "IS extraction and post-processing in progress...\n\n";

$genomeVector = $config->get_value("genomeVector");
my $vectorString = $config->get_value("vectorString");
my $suffOut = $config->get_value("suffOut");
if ($type eq "AGILENT") {
	if ($paired){
		my $genomeVectorIndexBlat = $config->get_value("genomeVectorIndexBlat");
		my $blatAligner = $config->get_value("blatAligner");
		my $minIden = $config->get_value("minIden");
		my $range = $config->get_value("range");
		my $samtools = $config->get_value("samtools");
		my $specificity = $config->get_value("specificity");
		system "(echo \"perl -I $libDir $scriptDir/extractIS.pl -l $libDir -p -aIn $alignmentOut.sorted.bam -o $output_dir -s $scriptDir -r $genomeVector -i $genomeVectorIndexBlat -v $vectorString -t $type   -bla $blatAligner -minIden $minIden -range $range -samtools $samtools -specificity $specificity \" )";
		system("perl -I $libDir $scriptDir/extractIS.pl -l $libDir -p -aIn $alignmentOut.sorted.bam -o $output_dir -s $scriptDir -r $genomeVector -i $genomeVectorIndexBlat -v $vectorString -t $type -bla $blatAligner  -minIden $minIden -range $range -samtools $samtools -specificity $specificity");
	
	} else {
	        # for Single End raeds.
	        my $vectorString  = $config->get_value("vectorString"); 
	        my $vectorIndexOut = $config->get_value("vectorIndexOut");
	        my $genomeIndexOut = $config->get_value("genomeIndexOut");
	        my $alignerBWA = $config->get_value("alignerBWA");
		my $minIden = $config->get_value("minIden");
	        my $range = $config->get_value("range");
	     	if ($config->has_value("additionalIndexSingleEnd")){
       	                my $genomeVectorIndexBlat = $config->get_value("genomeVectorIndexBlat");
			my $blatAligner = $config->get_value("blatAligner");
		        system("perl -I $libDir $scriptDir/extractIS.pl -aIn $workingDir/$alignmentOut.sam -o $workingDir -s $scriptDir -v $vectorString  -vecFile $workingDir/$vectorIndexOut  -t $type -genFile   $workingDir/$genomeIndexOut -fOut $workingDir/$forwardOut.fastq -aBWA $alignerBWA   -i $genomeVectorIndexBlat   -bla $blatAligner  ");
	        } else {
                	my $genomeVectorIndexBlat = $config->get_value("genomeVectorIndexBlat");
			my $blatAligner = $config->get_value("blatAligner");
			system("echo \"perl -I $libDir $scriptDir/extractIS.pl -aIn $alignmentOut.sam -o $output_dir -s $scriptDir -v $vectorString  -vecFile $vectorIndexOut  -genFile   $genomeIndexOut -fOut $suffOut.fastq -aBWA $aligner -t $type  -i $genomeVectorIndexBlat -bla $blatAligner -minIden $minIden -range $range\"  ");
			system("perl $scriptDir/extractIS.pl -aIn $alignmentOut.sam -o $output_dir -s $scriptDir -v $vectorString  -vecFile $vectorIndexOut  -genFile   $genomeIndexOut -fOut $suffOut.fastq -aBWA $aligner -t $type  -i $genomeVectorIndexBlat -bla $blatAligner  -minIden $minIden -range $range ");
	        }

	}
} else {
	my $alScore = 0.9;
	my $idClus = 0.9;
	my $range = 10;
	my $notMatchThreshold= 1000;
	my $minIden= 95;
	
	if ($config->has_value("alScore")){
		$alScore = $config->get_value("alScore");
	}
	if ($config->has_value("idClus")) {
		$idClus = $config->get_value("idClus");
	}
	if ($config->has_value("range")) {
		$range = $config->get_value("range");
	}

        if ($config->has_value("notMatchThreshold")) {
                $notMatchThreshold = $config->get_value("notMatchThreshold");
        }
        if ($config->has_value("minIden")) {
                $minIden = $config->get_value("minIden");
        }

	my $genomeVectorIndexBlat = $config->get_value("genomeVectorIndexBlat");
	my $blatAligner = $config->get_value("blatAligner");

	if (($lamType eq "both" || $lamType eq "extracting")) {
		system("echo \"perl $scriptDir/extractIS.pl  -o $output_dir -s $scriptDir -t $type -m $lenLam -a $alScore -idClus $idClus -range $range -notMatchThreshold  $notMatchThreshold -bla  $blatAligner   -i  $genomeVectorIndexBlat -minIden $minIden\" ");
		system("perl $scriptDir/extractIS.pl  -o $output_dir -s $scriptDir -t $type -m $lenLam -a $alScore -idClus $idClus -range $range  -notMatchThreshold  $notMatchThreshold    -bla   $blatAligner   -i $genomeVectorIndexBlat -minIden $minIden ");
	}
	else {
		print "No extraction step required\n"
	}
}

##############################################################################################
###############################################################################################
#Approx IS extraction (optional)
if ($type eq "AGILENT") {
        if ($paired){
                my $approxIS = $config->get_value("approxIS");
                        if ($approxIS eq "TRUE") {
                                print "##############################################################\n";
                                print "Approx IS extraction ...\n\n";
                                my $vectorString = $config->get_value("vectorString");
                                my $alignerBWA = $config->get_value("aligner");
                                my $refGenIndexBWA = $config->get_value("genomeIndexOut");
                                my $vecIndexBWA = $config->get_value("vectorIndexOut");
				
                                system ("bash $scriptDir/approxISforTS.sh  $output_dir $scriptDir $output_dir/completAlignment.sam $output_dir/completAlignment.nodup.bam $vectorString $refGenIndexBWA $vecIndexBWA $alignerBWA $samtools  ");
				my $UCSCAnnoFile = $config->get_value("UCSCAnnoFile");
				my $bedTools = $config->get_value("bedTools");
			#	print "IS annotation TS approx...\n\n";	
                                system("perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1  $output_dir/approx.resultsNoDup.csv.total.results ");
                              	system("bash $scriptDir/approxFormatTS.sh $output_dir  ");
				};
                };
};
###############################################################################################
#Extra filtering optional
if ($type eq "AGILENT") {
	my $extraFilt = $config->get_value("extraFilt");
		if ($extraFilt eq "TRUE") {
			my $genomeVectorIndexBlat = $config->get_value("genomeVectorIndexBlat");
			my $blatAligner = $config->get_value("blatAligner");
			print "###################################################################\n";
			print "Extra filtering in process...\n\n";
			system("echo bash  $scriptDir/extraFilteringTS.sh $vectorString   $genomeVectorIndexBlat $output_dir  $scriptDir $blatAligner ");
		 	system("bash  $scriptDir/extraFilteringTS.sh $vectorString   $genomeVectorIndexBlat $output_dir  $scriptDir $blatAligner ");
	
			print "###################################################################\n";
		#	print "Multiple Aligned reads processing...\n\n";
			my $alScore = $config->get_value("alScore");
		#	my $blatAligner = $config->get_value("blatAligner");
		#	my $genomeVectorIndexBlat = $config->get_value("genomeVectorIndexBlat");
			my $minIden = $config->get_value("minIden");
			my $range = $config->get_value("range");
			print "Multiple aligned reads processing TS ...\n\n";
			system("echo bash  $scriptDir/repeatsExtractTS.sh   $vectorString  $output_dir  $scriptDir $alScore $blatAligner $genomeVectorIndexBlat $minIden $range ");
			system("bash  $scriptDir/repeatsExtractTS.sh   $vectorString  $output_dir  $scriptDir $alScore $blatAligner $genomeVectorIndexBlat $minIden $range ");
	

			print "###################################################################\n";
			print "IS annotation TS...\n\n";
       			my $UCSCAnnoFile = $config->get_value("UCSCAnnoFile");
     			my $bedTools = $config->get_value("bedTools");
		        system("echo \"perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1 $output_dir/resultsNoDup.csv.total.results\" ");

			system("echo \"perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1 $output_dir/repeats.resultsNoDup.csv.total.results\" ");

			system("perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1  $output_dir/repeats.resultsNoDup.csv.total.results ");
     			system ("bash $scriptDir/repeatsFormatTS.sh $output_dir  ");


        		system("perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1  $output_dir/resultsNoDup.csv.total.results ");
	
	
}


		if ($extraFilt eq "FALSE") {
			print "\nNo extra filtering...\n";
			print "###################################################################\n";
       		 	print "Multiple aligned reads processing TS ...\n\n";
			my $alScore = $config->get_value("alScore");
                        my $blatAligner = $config->get_value("blatAligner");
                        my $genomeVectorIndexBlat = $config->get_value("genomeVectorIndexBlat");
                        my $minIden = $config->get_value("minIden");
                        my $range = $config->get_value("range");
        		system("echo bash  $scriptDir/repeatsExtractTS.sh   $vectorString  $output_dir  $scriptDir $alScore $blatAligner $genomeVectorIndexBlat $minIden $range");
     		   	system("bash  $scriptDir/repeatsExtractTS.sh   $vectorString  $output_dir  $scriptDir $alScore $blatAligner $genomeVectorIndexBlat $minIden $range ");
			print "###################################################################\n";
			print "IS annotation TS...\n\n";
			my $UCSCAnnoFile = $config->get_value("UCSCAnnoFile");
			my $bedTools = $config->get_value("bedTools");
			system("echo \"perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1 $output_dir/resultsNoDup.csv.total.results\" ");
			system("echo \"perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1 $output_dir/repeats.resultsNoDup.csv.total.results\" ");

			system("perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1  $output_dir/repeats.resultsNoDup.csv.total.results ");
			system ("bash $scriptDir/repeatsFormatTS.sh $output_dir  ");

			system("perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1  $output_dir/resultsNoDup.csv.total.results " );
}
}


#FOR LAM_PCR
print "###################################################################\n";
if ($config->has_value("lamType")) {
        my $alScore = 0.9;
        my $idClus = 0.9;
        my $range = 10;

        if ($config->has_value("alScore")){
               $alScore = $config->get_value("alScore");
       }

        if ($config->has_value("idClus")) {
               $idClus = $config->get_value("idClus");
        }
        if ($config->has_value("range")) {
               $range = $config->get_value("range");
        }

	if (($lamType eq "both" || $lamType eq "extracting")) {
        print "\nMultiple aligned reads processing LAM ...\n\n";
		system("echo \"perl $scriptDir/repeatsExtractLAM.extractIS.pl  -o $output_dir -s $scriptDir -t $type -m $lenLam -a $alScore -idClus $idClus -range $range\" ");
		system("perl $scriptDir/repeatsExtractLAM.extractIS.pl  -o $output_dir -s $scriptDir -t $type -m $lenLam -a $alScore -idClus $idClus -range $range ");	

		print "###################################################################\n";
    		print "IS annotation LAM...\n\n";
	        my $UCSCAnnoFile = $config->get_value("UCSCAnnoFile");
	        my $bedTools = $config->get_value("bedTools");
        	system("echo \"perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1 $output_dir/resultsNoDup.csv.total.results.r\" ");
		
	        system("echo \"perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1 $output_dir/resultsNoDup.csv.total.results\" ");
         	system("perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1  $output_dir/resultsNoDup.csv.total.results.r ");
	        system ("bash $scriptDir/repeatsFormatTS.sh $output_dir ");

        	system("perl -I $libDir $scriptDir/annotation.pl -o $output_dir -s $scriptDir -t $bedTools -a1 $UCSCAnnoFile  -r1  $output_dir/resultsNoDup.csv.total.results ");

}
}

###############################################################################################
###############################################################################################
print "###################################################################\n";

if ($type eq "LAM"){
	if (($lamType eq "both" || $lamType eq "extracting")) {	
		my $idClus = 0.9;
		print "Generating General Statistics ...\n\n";
		my $runName = $config->get_value("runName");

		my $len = $config->get_value("len");
		my $alScore = $config->get_value("alScore");
	        if ($config->has_value("idClus")) {
               $idClus = $config->get_value("idClus");
        }

#		my $idClus = $config->get_value("idClus");
		my $range = $config->get_value("range");
		my $notMatchThreshold = $config->get_value("notMatchThreshold");

		system ("bash $scriptDir/generalStatisticsLAM.sh $output_dir $forward $len  $alScore  $idClus  $range   $notMatchThreshold ");
		system ("bash $scriptDir/finalResultsFormatLAM.sh $output_dir $runName $len  $alScore  $idClus  $range   $notMatchThreshold ");
}

}

if ($type eq "AGILENT") {
print "Generating General Statistics ...\n\n";
	 if ($paired){
		 my $approxIS = $config->get_value("approxIS");
                        if ($approxIS eq "FALSE") {

				my $sampleName = $config->get_value("sampleName");
				system ("bash $scriptDir/generalStatisticsSureSelectPair.sh $output_dir $forward $vectorString $samtools ");
				system ("bash $scriptDir/finalResultsFormatTS.sh $output_dir $sampleName $samtools");
			}
			
			if ($approxIS eq "TRUE") {

                                my $sampleName = $config->get_value("sampleName");
                                system ("bash $scriptDir/approx.generalStatisticsSureSelectPair.sh $output_dir $forward $vectorString $samtools");
                                system ("bash $scriptDir/approx.finalResultsFormatTS.sh $output_dir $sampleName ")
                        }
	
	} else {
		my $sampleName = $config->get_value("sampleName");
		system ("bash $scriptDir/generalStatisticsSureSelectSingle.sh  $output_dir $forward $vectorString $samtools ");
		system ("bash $scriptDir/finalResultsFormatTS.sh $output_dir $sampleName $samtools ")
	}
}	
print "Finished.\n";
print "###################################################################\n";
###############################################################################################
###############################################################################################
###############################################################################################	
