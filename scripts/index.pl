#!/usr/bin/perl -w

use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use GeneisLib;
use Scalar::Util qw(looks_like_number);
use Scalar::Util qw(isvstring);


my @usage;
push @usage, "\nUsage:  index.pl  <-g genome>  <required>  <-v virus>  <required> <-a aligner>  <required> [options] \n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this information.\n";
push @usage, "  -g, --genome  Genome fasta file <required>.\n";
push @usage, "  -v, --vector   Virus fasta file <required>.\n";
push @usage, "  -a, --aligner   Aligner full path, bwa/bowtie2 <required>.\n";
push @usage, "  -iOut, --indexOut   Name for output indexed files.Default is refVec\n";

push @usage, "  -o, --output   Full path of a directory to store results, default is current working directory.\n\n";

push @usage, "  This program provides indexing of reference genome:\n";
push @usage, "#########################################################################\n\n";


###############################################################################################

my $help;
my $g_file;
my $v_file;

my $a_file;
my $iOut_value;
my $output_dir;


GetOptions
(
 'h|help|?'    => \$help,
 'g=s'    => \$g_file,
 'v=s'    => \$v_file,
 'a=s'    => \$a_file,
 'iOut=s'    => \$iOut_value,
 'output=s'    => \$output_dir,
);

###############################################################################################
if ($help) {
   print @usage;
   exit(0);
}

if (defined $g_file) {
   if (!-e "$g_file.fa"){
   	print "\nThe genome fasta file $g_file.fa does not exist!\n\n";
        print @usage;
	exit;
   }
}else{
    print "Please provide a genome fasta file!\n\n";
    print @usage;
    exit;
}



if (defined $v_file) {
   if (!-e "$v_file.fa"){
   	print "\nThe virus fasta file $v_file.fa does not exist!\n\n";
        print @usage;
	exit;
   }
}else{
    print "Please provide a virus fasta file!\n\n";
    print @usage;
    exit;
}



#if (defined $a_file) {
#   if (!-e $a_file){
#   	print "\nThe aligner $a_file does not exist!\n\n";
#        print @usage;
#	exit;
#   }
#}else{
#    print "Please provide a right aligner path!\n\n";
#    print @usage;
#    exit;
#}


if (defined $iOut_value) {
    if (isvstring $iOut_value){
    	print "\nThe $iOut_value is name for output indexed files.\n";
        
    }
}else{
    $iOut_value = refVec;
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

###############################################################################################
## Indexing

print "\n Indexing in progress...\n" ;
print "\nResults will be stored in  $output_dir\n\n";

# Construct the reference fasta files for all the aligner

system ("cat  $g_file.fa $v_file.fa > $output_dir/$iOut_value.fa");

if ($a_file =~ /bwa/) {
      system("$a_file index -a bwtsw $output_dir/$iOut_value.fa");

}
elsif ( $a_file=~ /blast/ ){
      system("formatdb -i $output_dir/$iOut_value.fa -p F -n $output_dir/$iOut_value");
}
else {
     #this is just to concatenate "build" word with bowtie2 as bowtie2-build is required in command.
     my $b_file = "$a_file-build";
     #system "echo $b_file";
     system("$b_file $output_dir/$iOut_value.fa $output_dir/$iOut_value");
}
##############################################################################################
