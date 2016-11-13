#!/usr/bin/perl -w


use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Scalar::Util qw(looks_like_number);
use Scalar::Util qw(isvstring);


my @usage;
push @usage, "\nUsage:  sortBarcode.pl  <-f forward file>  <-r reverse file> <-tt tag table> <-o results folder> [options] \n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this infrOutmation.\n";
push @usage, "  -f, --forward   Forward FASTQ file <required>.\n";
push @usage, "  -r, --reverse   Reverse FASTQ file.\n";
push @usage, "  -tt, --tagTable   Full path of cutadapt tool is needed <required>. \n";
push @usage, "  -g, --group   Number of reads in concurrent thread. Default is 5000000\n";
push @usage, "  -l, --len   minimal fragment length after sorting. Default is 25\n";
push @usage, "  -l, --len   minimal quality of read. Default is 30\n";
push @usage, "  -p, --proc   Concurrent number of processing threads. Default is 1\n";
push @usage, "  -lnk, --linker   Adapter for forward file.Default is CCTAACTGCTGTGCCACT\n";
push @usage, "  -mega, --megaprimer   Adapter for reverse file.Default is GATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA\n";
push @usage, "  -sOut, --suffOut   Name for suffix output file. Default is filtTrim\n";
push @usage, "  -anchor --anchorString   Sequence after megaprimer used to verify the trueness of the read. len(anchor)<30\n";
push @usage, "  -anchorMM --anchorStringMM   Mismatches allowed in anchor. Default is 2.\n";
push @usage, "  -s, --scripts   Full path of a directory that store scripts files.\n\n";
push @usage, "  -o, --output   Full path of a directory to store results.Default is current working directory.\n\n";
push @usage, "  -readHeader --readHeader   Raed Header\n";
push @usage, "This program sort a duoble barcoded set of reads with the provided FASTQ files and tagTable\n\n";
push @usage, "#########################################################################\n\n";



###############################################################################################

my $help;
my $f_file;
my $r_file;
my $tagTable;
my $len;
my $qual;
my $linker;
my $megaPrimer;
my $Out_value;
my $group;
my $scripts_dir;
my $output_dir;
my $threads;
my $sorting;
my $anchor;
my $anchorMM;
my $readHeader;

GetOptions
(
 'h|help|?'    => \$help,
 'f=s'    => \$f_file,
 'r=s'    => \$r_file,
 'tt=s'    => \$tagTable,
 'g=s'    => \$group,
 'lnk=s'    => \$linker,
 'mega=s'    => \$megaPrimer,
 'l=s'    => \$len,
 'qual=s'    => \$qual,
 'p=s'    => \$threads,
 'sOut=s'    => \$Out_value,
 'output=s'    => \$output_dir,
 's=s'    => \$scripts_dir,
 'ty=s'	=> \$sorting,
 'anchor=s'  => \$anchor,
 'anchorMM=s'  => \$anchorMM,
 'readHeader=s'    => \$readHeader,
 'readHeader=s'    => \$readHeader,
 'sk=s'    => \$skewer,
);

print "checking $threads";
print "check $readHeader";
###############################################################################################
if ($help) {
   print @usage;
   exit(0);
}

if (defined $f_file) {
   if (!-e $f_file){
   	print "\nThe forward file $f_file does not exist!\n\n";
        print @usage;
	exit;
   }
}else{
    print "Please provide a forward file!\n\n";
    print @usage;
    exit;
}



if (defined $r_file) {
   if (!-e $r_file){
   	print "\nThe reverse file $r_file does not exist!\n\n";
        print @usage;
	exit;
   }
}


if (defined $tagTable) {
   if (!-e $tagTable){
   	print "\nThe tag table $tagTable does not exist.\n";
	print @usage;
	exit;
   }
}



if (defined $len) {
   if (looks_like_number($len)){
   	print "\nMinimal read length is $len.\n";
   }
}else{
    $len = 25;
}

if (defined $qual) {
   if (looks_like_number($qual)){
        print "\nMinimal read length is $qual.\n";
   }
}else{
    $qual = 30;
}

if (defined $threads) {
   if (looks_like_number($threads)){
   	print "\nConcurrent threads = $threads.\n";
   }
}else{
    $threads = 1;
}


if (defined $readHeader) {
    if (isvstring ($readHeader)){
        print "\nThe $readHeader is being used for forward file.\n";
        print @usage;
        exit;
    }
}

if (defined $linker) {
    if (isvstring ($linker)){
    	print "\nThe $linker is being used for forward file.\n";
        print @usage;
	exit;
    }
}


if (defined $sorting) {
	if ($sorting ne "both" &&  $sorting ne "sorting" && $sorting ne "extracting") {
		$sorting="both";
	}
} else {
	$sorting="both";
}


if (defined $megaPrimer) {
    if (isvstring ($megaPrimer)){
    	print "\nThe $megaPrimer is being used for reverse file.\n";
        print @usage;
	exit;
    }
}



if (defined $Out_value) {
    if (isvstring $Out_value){
    	print "\n$Out_value is the prefix ot the LAM files.\n";
        
    }
}else{
    $Out_value = MiS;
}


if (defined $output_dir) {
    if (!-e $output_dir){
    	print "\nThe output $output_dir does not exist!\n\n";
        print @usage;
	exit;
    }
}else{
    $output_dir = getcwd;
}

if (defined $scripts_dir) {
    if (!-d $scripts_dir){
    	print "\nThe folder $scripts_dir does not exist!\n\n";
        print @usage;
	exit;
    }
}else{
    	print "\nThe script folder must be provided!\n\n";
        print @usage;
	exit;
}

if (!-d "$output_dir/tmpR") {
	system "(mkdir $output_dir/tmpR)";
    	print "\nThe temporary directory $output_dir/tmpR is created\n\n";
}

if (!-d "$output_dir/tmpF") {
	system "(mkdir $output_dir/tmpF)";
    	print "\nThe temporary directory $output_dir/tmpF is created\n\n";
}


###############################################################################################
# Quality filtering and adapter trimmming.
print "$readHeader";
print "\nResults will be stored in  $output_dir\n";

if ($sorting eq "both" || $sorting eq "sorting") {
	print "\nForward Splitting in progress...\n" ;
	system "(echo \"zcat $f_file| split -l $group - $output_dir/tmpF/tmp\")";
	system "(       zcat $f_file| split -l $group - $output_dir/tmpF/tmp)";
	print "\nReverse Splitting in progress...\n" ;
	system "(echo \"zcat $r_file| split -l $group - $output_dir/tmpR/tmp\")";
	system "(       zcat $r_file| split -l $group - $output_dir/tmpR/tmp)";
	print "\nStop Splitting\n";
	system ("echo \"bash $scripts_dir/sorting.sh $tagTable $output_dir $scripts_dir $threads $readHeader\"");
	system ("       bash $scripts_dir/sorting.sh $tagTable $output_dir $scripts_dir $threads $readHeader");
	print "\nEnd Sorting\nStarted Merging/Trimming\n";
}
	
if ($sorting eq "both" || $sorting eq "extracting") {
#if ($sorting eq "both" || $sorting eq "sorting") {
	print "$anchor\n";
	if ($anchor eq "noAnchor") {
		system("echo \"bash $scripts_dir/merging.sh $output_dir $Out_value $qual $skewer \"");
		system("       bash $scripts_dir/merging.sh $output_dir $Out_value $qual $skewer ");
	}
	else {
		if (!defined $anchorMM) {
			$anchorMM=0.1
		}
		system("echo \"bash $scripts_dir/merging.sh $output_dir $Out_value $qual $skewer $anchor $anchorMM\"");
		system("       bash $scripts_dir/merging.sh $output_dir $Out_value $qual $skewer $anchor $anchorMM");
	}
	#system("echo \"bash $scripts_dir/merging.sh $megaPrimer $linker $output_dir $Out_value\"");
	#system("       bash $scripts_dir/merging.sh $megaPrimer $linker $output_dir $Out_value");
	print "End Merging/Trimming\n";
}

################################################################################################
