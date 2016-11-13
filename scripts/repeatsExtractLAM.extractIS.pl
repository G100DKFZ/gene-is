#!/usr/bin/perl -w

use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Scalar::Util qw(looks_like_number);
use Scalar::Util qw(isvstring);


my @usage;
push @usage, "\nUsage:  LAM repeats extractIS.pl  [options] \n\n";
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
push @usage, "  -fOut, --forward filtered and trimmed FASTQ file. (REQUIRED ONLY FOR SINGLE END READS)\n\n";
push @usage, "  -aBWA, --path to BWA aligner is needed.(REQUIRED ONLY FOR SINGLE END READS)\n\n";
push @usage, "  -idClus, --minimum identity percentage within reads clusters in LAM (0.9)\n\n";
push @usage, "  -range, --maximum topographical distance to collapse IS (10 bp)\n\n";

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
GetOptions
(
 'h|help|?'    => \$help,
# 'p'    => \$paired,
# 'aIn=s'    => \$inAlnFile,
 'o=s'    => \$output_dir,
 's=s'    => \$scriptDir,
 'a=s'    => \$alScore,
 'l=s'    => \$libDir,
 'm=s'    => \$len,
 #'t=s'    => \$type,
 #'r=s'    => \$refGenome,
 #'i=s'    => \$refGenomeIndexBlat,
# 'v=s'    => \$vectorString,
# 'vecFile=s'    => \$vectorFile,
# 'genFile=s'    => \$genomeFile,
# 'fOut=s'    => \$forwardOut,
# 'aBWA=s'    => \$alignerBWA,
 'idClus=s'  => \$idClus,
 'range=s'  => \$range,
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

if (!defined $idClus) {
	$idClus=0.9
}

if (!defined $range) {
	$range=10
}

###############################################################################################
# IS extraction for MH

#print "\nIS extraction for LAM repeats/multipleAligned in progress...\n" ;
#print "\nResults will be stored in  $output_dir\n";

if (defined $alScore) {
	system ("echo \"bash $scriptDir/repeatsExtractLAM.extractLAMIS.sh $output_dir $scriptDir $len $alScore $idClus $range\""); 
	system ("bash $scriptDir/repeatsExtractLAM.extractLAMIS.sh $output_dir $scriptDir $len $alScore $idClus $range") ;
} else {
	system ("echo \"bash $scriptDir/repeatsExtractLAM.extractLAMIS.sh $output_dir $scriptDir $len 0 $idClus $range\""); 
	system ("bash $scriptDir/repeatsExtractLAM.extractLAMIS.sh $output_dir $scriptDir $len 0 $idClus $range") ;
	}


################################################################################################
