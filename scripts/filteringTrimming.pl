#!/usr/bin/perl -w

use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Scalar::Util qw(looks_like_number);
use Scalar::Util qw(isvstring);

my @usage;
push @usage, "\nUsage:  filterTrimming.pl  <-f forward file>  <required>  <-skewer full path is required> <required>. [options] \n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this infrOutmation.\n";
push @usage, "  -f, --forward   Forward FASTQ file <required>.\n";
#push @usage, "  -cutadapt, --cutadaptTrimmer   Full path of cutadapt tool is needed <required>. \n";
push @usage, "  -sk, --skewer   Full path of skewer tool is needed <required>. \n";
push @usage, "  -r, --reverse   Reverse FASTQ file.\n";
push @usage, "  -qual, --quality   Quality value for filteration <0-40>.Default is 20\n";
#push @usage, "  -perc, --percentage Percentage value for filteration <0-100>. Default is 80\n";
push @usage, "  -adaptF, --adapterForward   Adapter for forward file.Default is GATCGGAAGAGCACACGTCTGAACTCCAGTCAC\n";
push @usage, "  -adaptR, --adapterReverse   Adapter for reverse file.Default is AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\n";
push @usage, "  -sOut, --suffOut   Name for suffix output file. Default is filtTrim\n";
push @usage, "  -o, --output   Full path of a directory to store results.Default is current working directory.\n\n";
push @usage, "This program quality filter and trim adapters of the provided FASTQ files\n\n";
push @usage, "#########################################################################\n\n";

###############################################################################################

my $help;
my $f_file;
my $r_file;
my $qual_value;
#my $perc_value;
my $adaptF_value;
my $adaptR_value;
my $Out_value;
my $cutadapt_file;
my $output_dir;

GetOptions
(
 'h|help|?'    => \$help,
 'f=s'    => \$f_file,
 'r=s'    => \$r_file,
 'qual=s'    => \$qual_value,
# 'perc=s'    => \$perc_value,
 'adaptF=s'    => \$adaptF_value,
 'adaptR=s'    => \$adaptR_value,
 'cutadapt=s'    => \$cutadapt_file,
 'sOut=s'    => \$Out_value,
 'output=s'    => \$output_dir,
 'sk=s'    => \$skewer_dir,
);


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

#if (defined $cutadapt_file) {
#   if (!-e $cutadapt_file){
#   	print "\nThe forward file $cutadapt_file does not exist!\n\n";
#        print @usage;
#	exit;
#   }
#}else{
#    print "Please provide a cutadapt tool path !\n\n";
#    print @usage;
#    exit;
#}

if (defined $r_file) {
   if (!-e $r_file){
   	print "\nThe reverse file $r_file does not exist!\n\n";
        print @usage;
	exit;
   }
}

if (defined $qual_value) {
   if (looks_like_number($qual_value)){
   	print "\nQuality value is $qual_value.\n"
   }
}else{
    $qual_value = 20;
}

#if (defined $perc_value) {
 #  if (looks_like_number($perc_value)){
  # 	print "\nQuality percentage value is $perc_value.\n"
  # }
#}else{
 #   $perc_value = 80;
#}
if (defined $adaptF_value) {
    if (isvstring ($adaptF_value)){
    	print "\nThe $adaptF_value is being used for forward file.\n";
        print @usage;
	exit;
    }
}else{
    $adaptF_value = GATCGGAAGAGCACACGTCTGAACTCCAGTCAC ;
}

if (defined $adaptR_value) {
    if (isvstring ($adaptR_value)){
    	print "\nThe $adaptR_value is being used for reverse file.\n";
        print @usage;
	exit;
    }
}else{
    $adaptR_value = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT ;
}

if (defined $Out_value) {
    if (isvstring $Out_value){
    	print "\nThe $Out_value is name for output forward file.\n";
        
    }
}else{
    $Out_value = filtTrim;
}

if (defined $output_dir) {
    if (!-e $output_dir){
    	print "\nThe output directrOuty $output_dir does not exist!\n\n";
        print @usage;
	exit;
    }
}else{
    $output_dir = getcwd;
}

###############################################################################################
#print "\nQuality filtering and adapter trimming in progress...\n" ;
print "\nResults will be stored in  $output_dir\n";

if (defined $r_file) {
   

	#system "(fastq_quality_fiilter  -q $qual_value -p $perc_value -i $f_file -o $output_dir/filt11.fastq   -Q33)";
	system "(echo \"$skewer_dir -x $adaptF_value -y $adaptR_value -q $qual_value -l 50 -o $output_dir/$Out_value $f_file $r_file\")";
	system "($skewer_dir -x $adaptF_value -y $adaptR_value -q $qual_value -l 20 -o $output_dir/$Out_value $f_file $r_file)";
	#system "(fastq_quality_filter  -q $qual_value -p $perc_value -i $r_file -o $output_dir/filt22.fastq  -Q33)";
}else{
	system "(echo \"$skewer_dir -x $adaptF_value -q $qual_value -l 50 -o $output_dir/$Out_value $f_file\")"; 
	system "($skewer_dir -x $adaptF_value -q $qual_value -l 20 -o $output_dir/$Out_value $f_file)";

}
################################################################################################

