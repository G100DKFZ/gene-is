#!/usr/bin/perl -w

use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Scalar::Util qw(looks_like_number);
use Scalar::Util qw(isvstring);


my @usage;
push @usage, "\nUsage:  alignment.pl  <-f forward file>  <required> <-gv indexed files>  <required> <-a aligner>  <required>   [options] \n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this information.\n";
push @usage, "  -f, --forward   Forward FASTQ file <required>.\n";
push @usage, "  -a, --aligner   Aligner full path, bwa/bowtie2 <required>.\n";
push @usage, "  -gv, --genomeVector Name of indexed reference file    <required>.\n";
push @usage, "  -r, --reverse   Reverse FASTQ file.\n";
push @usage, "  -aOut, --alignmentOut  Name for output sam file.Default is completAlignment\n";
push @usage, "  -p, --parallel   Set number of parallel threads.\n";
push @usage, "  -t, --type   Analysis type (AGILENT|LAM). Dafault [LAM].\n";
push @usage, "  -o, --output   Full path of a directory to store results, default is current working directory.\n";
push @usage, "  -s, --scripts   Full path of a directory where are the scripts.\n";
push @usage, "  This program provides the alignment bam file:\n";
push @usage, "#########################################################################\n\n";

###############################################################################################

my $help;
my $f_file;
my $r_file;
my $a_file;
my $paired;
my $gv_file;
my $aOut_file;
my $output_dir;
my $scripts_dir;
my $threads;
my $type;
my $samtools;
GetOptions
(
 'h|help|?'    => \$help,
 'f=s'    => \$f_file,
 'r=s'    => \$r_file,
 'a=s'    => \$a_file,
 'gv=s'    => \$gv_file,
 'aOut=s'    => \$aOut_value,
 'output=s'    => \$output_dir,
 's=s'    => \$scripts_dir,
 'p=s'		=> \$threads,
 't=s'		=> \$type,
 'sam=s'          => \$samtools,
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

if ($type eq "AGILENT") {
	if (defined $f_file) {
   		if (!-e "$output_dir/$f_file"){
		   	print "\nThe forward file $output_dir/$f_file does not exist!\n\n";
		        print @usage;
			exit;
		   }
	}else{
	    print "Please provide a forward file!\n\n";
	    print @usage;
	    exit;
	}
}

if (defined $r_file) {
   $paired=1;
   if (!-e "$output_dir/$r_file"){
   	print "\nThe reverse file $output_dir/$r_file does not exist!\n\n";
        print @usage;
	exit;
   }
}


if (defined $a_file) {
   if (!-e $a_file){
   	print "\nThe aligner $a_file does not exist!\n\n";
        print @usage;
	exit;
   }
}else{
    print "Please provide a right aligner path!\n\n";
    print @usage;
    exit;
}


if (defined $gv_file) {
   if ((!-e "$gv_file") && (!-e "$gv_file.fa")){
   	print "\nThe indexed file $gv_file does not exist!\n\n";
        print @usage;
	exit;
   }
}else{
    print "Please provide a right indexed file!\n\n";
    print @usage;
    exit;
}

if (defined $aOut_value) {
    if (isvstring $aOut_value){
    	print "\nThe $aOut_value is name for alignment sam file.\n";
        
    }
}else{
    $aOut_value = completeAlignment;
}  

###############################################################################################
#print "\n Alignment in progress...\n" ;
#print "\nResults will be stored in  $output_dir\n";
print "Alignemnt type $type";

####### TS Section ###########
if ($type eq "AGILENT") {
# For BWA aligner
	if ($a_file =~ /bwa/) {
		print "\nBWA is used as Aligner\n\n";
		if (defined $r_file) {
       	       	# For paired end
       	 		system("echo \"$a_file mem -M -t $threads $gv_file  $output_dir/$f_file  $output_dir/$r_file  > $output_dir/$aOut_value.sam \" ");
       	 		system("$a_file mem -M -t $threads $gv_file  $output_dir/$f_file  $output_dir/$r_file  > $output_dir/$aOut_value.sam ");
			#system("echo \"awk  -F $scripts_dir/correctSam.awk $output_dir/$aOut_value.sam > app.sam;mv app.sam  $output_dir/$aOut_value.sam\"");
			#`awk -f $scripts_dir/correctSam.awk $output_dir/$aOut_value.sam > app.sam;mv app.sam  $output_dir/$aOut_value.sam`;
			`$samtools view -Sb $output_dir/$aOut_value.sam  > $output_dir/$aOut_value.bam `;
			`$samtools sort $output_dir/$aOut_value.bam $output_dir/$aOut_value.sorted `;
			`$samtools index $output_dir/$aOut_value.sorted.bam `;
			`$samtools rmdup $output_dir/$aOut_value.sorted.bam $output_dir/$aOut_value.nodup.bam `;
			`$samtools index $output_dir/$aOut_value.nodup.bam `;
    		}else{
       	       		# For single end
        	     	system("echo \"$a_file mem -M  -t $threads $gv_file  $output_dir/$f_file  > $output_dir/$aOut_value.sam\" ");
      	       		system("$a_file mem -M  -t $threads $gv_file  $output_dir/$f_file  > $output_dir/$aOut_value.sam ");
			`$samtools view -Sb $output_dir/$aOut_value.sam > $output_dir/$aOut_value.bam `;
		       #`samtools sort $output_dir/$aOut_value.bam $output_dir/$aOut_value.sorted`;
                       # `samtools index $output_dir/$aOut_value.sorted.bam`;
                       # `samtools rmdup -s $output_dir/$aOut_value.sorted.bam $output_dir/$aOut_value.nodup.bam`;
                        `$samtools index $output_dir/$aOut_value.nodup.bam `;
       
            	}

# For Bowtie2 aligner
	}else{
		if ( $paired ) {
              	# For paired end
        		system("echo \"$a_file -p $threads --no-mixed --no-discordant --very-sensitive-local -x $gv_file  -1 $output_dir/$f_file -2 $output_dir/$r_file  -S $output_dir/$aOut_value.sam\"");
        		system("$a_file -p $threads --no-mixed --no-discordant --very-sensitive-local -x $gv_file  -1 $output_dir/$f_file -2 $output_dir/$r_file  -S $output_dir/$aOut_value.sam ");
    		}else{
              	# For single end
        		system("$a_file -p $threads --very-sensitive-local -x $gv_file  -q $output_dir/$f_file  -S $output_dir/$aOut_value.sam ");
      	 
        	}
		`$samtools view -Sb $output_dir/$aOut_value.sam > $output_dir/$aOut_value.bam`;
		`$samtools sort $output_dir/$aOut_value.bam $output_dir/$aOut_value.sorted`;
		`$samtools index $output_dir/$aOut_value.sorted.bam`;
		if ($paired) {
			`$samtools rmdup $output_dir/$aOut_value.sorted.bam $output_dir/$aOut_value.nodup.bam`;
		} else {
			`$samtools rmdup -s $output_dir/$aOut_value.sorted.bam $output_dir/$aOut_value.nodup.bam`;
		}		
	}
}


#################################
####### LAM Section ############
else {

	`echo "bash $scripts_dir/alignLAM.sh $output_dir $gv_file $a_file $threads $samtools"`;
	`bash $scripts_dir/alignLAM.sh $output_dir $gv_file $a_file $threads $samtools`;
}
###############################################################################################
