#!/usr/bin/perl -w

use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Scalar::Util qw(looks_like_number);
use Scalar::Util qw(isvstring);


my @usage;
push @usage, "\nUsage:  extractIS.pl  [options] \n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this infrOutmation.\n";
push @usage, "  -aIn, --alignmentIn   The file that contains alignment data. Should be in sorted BAM format in case of paired=true. In case of single end data file should be in SAM format\n";
push @usage, "  -o, --output_dir   Working directory. All working files should be placed in this directory.\n";

push @usage, "  -s, --scriptDir   Script directory.\n";
push @usage, "  -a, --alScore   Homology Score.\n";
push @usage, "  -l, --libDir   Perl Module directory.\n";
push @usage, "  -m, --len   reads minimum length.\n";
push @usage, "  -v, --vectorString   ID vector name.\n\n";
push @usage, "#########################################################################\n\n";
push @usage, "  -p, --paired   The data are from paired reads.\n";
push @usage, "  -r, --refGen   Reference genome.(REQUIRED ONLY FOR PAIRED END READS) \n";
push @usage, "  -bla, --blatAligner   path to BLAT aligner.\n";
push @usage, "#########################################################################\n\n";
push @usage, "  -vecFile, --vectorIndexedFiles (REQUIRED ONLY FOR SINGLE END READS).\n\n";
push @usage, "  -genFile, --genome indexed files (just genome index files without vector) (REQUIRED ONLY FOR SINGLE END READS) .\n";
push @usage, "  -fOut, --forward filtered and trimmed FASTQ file. (REQUIRED ONLY FOR SINGLE END READS)\n\n";
push @usage, "  -aBWA, --path to BWA aligner is needed.(REQUIRED ONLY FOR SINGLE END READS)\n\n";
push @usage, "  -idClus, --minimum identity percentage within reads clusters in LAM (0.9)\n\n";
push @usage, "  -range, --maximum topographical distance to collapse IS (10 bp)\n\n";
push @usage, "	-notMatchThreshold,  --maximum number of not match bases allowed before genomic match part\n\n";
push @usage, "  -minIden, --minimum identity percentage for re-alignment with LAM\n\n"; 
push @usage, "This program parse the alignement file and extract IS\n\n";
push @usage, "#########################################################################\n\n";



###############################################################################################

my $help;
my $paired;
my $len;
my $inAlnFile;
my $output_dir; 
my $scriptDir; 
my $libDir;
my $type;
my $alScore;
my $refGenome; 
my $refGenomeIndexBlat;
my $vectorString;
my $vectorFile;
my $genomeFile;
my $forwardOut;
my $alignerBWA;
my $idClus;
my $range;
my $notMatchThreshold;
my $minIden;
my $bla;
my $samtools;
my $specificity;
GetOptions
(
 'h|help|?'    => \$help,
 'p'    => \$paired,
 'aIn=s'    => \$inAlnFile,
 'o=s'    => \$output_dir,
 's=s'    => \$scriptDir,
 'a=s'    => \$alScore,
 'l=s'    => \$libDir,
 'm=s'    => \$len,
 't=s'    => \$type,
 'r=s'    => \$refGenome,
 'i=s'    => \$refGenomeIndexBlat,
 'v=s'    => \$vectorString,
 'vecFile=s'    => \$vectorFile,
 'genFile=s'    => \$genomeFile,
 'fOut=s'    => \$forwardOut,
 'aBWA=s'    => \$alignerBWA,
 'idClus=s'  => \$idClus,
 'range=s'  => \$range,
 'notMatchThreshold=s' =>  \$notMatchThreshold,
 'minIden=s' =>  \$minIden,
 'bla=s'    => \$blatAligner,
 'samtools=s'    => \$samtools,
 'specificity=s' =>  \$specificity,
);


###############################################################################################
if ($help) {
   print @usage;
   exit(0);
}

if (defined $output_dir) {
	if (!-d $output_dir) {
		print "\nThe working directory $output_dir does not exist!\n\n";
		print @usage;
		exit;
	}
}else{
	print "Please provide a working directory !\n\n";
	print @usage;
	exit;
}


if (defined $scriptDir) {
	if (!-d $scriptDir) {
		print "\nThe script directory $scriptDir does not exist!\n\n";
		print @usage;
		exit;
	}
}else{
	print "Please provide a script directory !\n\n";
	print @usage;
	exit;
}

if ($type eq "AGILENT") {
	if (defined $inAlnFile) {
	   if (!-e "$output_dir/$inAlnFile"){
	   	print "\nThe file $inAlnFile does not exist!\n\n";
	        print @usage;
		exit;
	   }
	}else{
	    print "Please provide a $inAlnFile  file!\n\n";
	    print @usage;
	    exit;
	}

	if (defined $vectorString) {
	   if (isvstring ($vectorString)){
	   	print "\nThe  $vectorString is provided $vectorString!\n\n";
	   }
	}else{
		print "Please provide a correct $vectorString !\n\n";
		print @usage;
		exit;
	}
}

if (!defined $idClus) {
	$idClus=0.9
}

if (!defined $range) {
	$range=10
}

if (!defined $notMatchThreshold) {
	$notMatchThreshold=1000
}

if (!defined $notMatchThreshold) {
        $notMatchThreshold=95
}

###############################################################################################
# IS extraction.

#print "\nIS extraction in progress...\n" ;
print "\nResults will be stored in  $output_dir\n";

if ($type eq "AGILENT") {

	##############################
	##Requirement for Paired End Only

	if ($paired){
		if (defined $refGenomeIndexBlat) {
			if (!-e "$refGenomeIndexBlat") {
	   			print "\nThe reference refGenomeIndexBlat $refGenomeIndexBlat does not exist!\n\n";
	        		print @usage;
				exit;
	   		}
		}
		if (defined $refGenome) {
	   		if (!-e "$refGenome"){
	   			print "\nThe reference refGenome $refGenome does not exist!\n\n";
	        		print @usage;
				exit;
	   		}
		}else{
	    		print "Please provide a reference genome !\n\n";
	    		print @usage;
	    		exit;
		}
		system("echo \"bash $scriptDir/extractIS.sh $output_dir $inAlnFile $scriptDir $refGenome $vectorString $blatAligner  $libDir $refGenomeIndexBlat $minIden $range $samtools $specificity \"");
		system("bash $scriptDir/extractIS.sh $output_dir $inAlnFile $scriptDir $refGenome $vectorString $blatAligner  $libDir $refGenomeIndexBlat $minIden $range $samtools $specificity ");
	
	##############################
	##Requirement for Single End Only
	
	} else {
                if (defined $refGenomeIndexBlat) {
                        if (!-e "$refGenomeIndexBlat") {
                                print "\nThe reference refGenomeIndexBlat $refGenomeIndexBlat does not exist!\n\n";
                                print @usage;
                                exit;
                        }
                }

		if (defined $vectorFile) {
	   		if (!-e $vectorFile){
	   			print "\nThe file $vectorFile does not exist!\n\n";
				print @usage;
				exit;
	   		}
		}
	
		if (defined $genomeFile) {
	   		if (!-e $genomeFile){
	   			print "\nThe file $genomeFile does not exist!\n\n";
				print @usage;
				exit;
	   		}
		}	
	
		if (defined $forwardOut) {
	   		if (!-e "$output_dir/$forwardOut"){
	   			print "\nThe file $output_dir/$forwardOut does not exist!\n\n";
				print @usage;
	                	exit;
			}	
		}
		system(" bash $scriptDir/extractSingleEndIS.sh $inAlnFile $output_dir $scriptDir $vectorString $vectorFile $genomeFile $forwardOut $blatAligner  $alignerBWA $refGenomeIndexBlat $minIden  $range ");
		system(" echo \"bash $scriptDir/extractSingleEndIS.sh $inAlnFile $output_dir $scriptDir $vectorString $vectorFile $genomeFile $forwardOut  $blatAligner  $alignerBWA $refGenomeIndexBlat $minIden  $range \"");
	}
} else {
	if (defined $alScore) {
		system ("echo \"bash $scriptDir/extractLAMIS.sh $output_dir $scriptDir $len $alScore $idClus $range $notMatchThreshold $blatAligner $refGenomeIndexBlat $minIden\""); 
		system ("bash $scriptDir/extractLAMIS.sh $output_dir $scriptDir $len $alScore $idClus $range $notMatchThreshold  $blatAligner $refGenomeIndexBlat $minIden") ;
	} else {
		system ("echo \"bash $scriptDir/extractLAMIS.sh $output_dir $scriptDir $len 0 $idClus $range $notMatchThreshold  $blatAligner $refGenomeIndexBlat $minIden\""); 
		system ("bash $scriptDir/extractLAMIS.sh $output_dir $scriptDir $len 0 $idClus $range $notMatchThreshold   $blatAligner $refGenomeIndexBlat $minIden") ;
	}
}

################################################################################################
