#!/usr/bin/perl -w

use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use Scalar::Util qw(looks_like_number);
use Scalar::Util qw(isvstring);

my @usage;
push @usage, "\nUsage:  annotation.pl  [options] \n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this information.\n";
push @usage, "  -o, --output_dir   Working directory. All working files should be placed in this directory.\n";
push @usage, "  -s, --scriptDir   Script directory.\n";
push @usage, "  -l, --libDir   Perl Module directory.\n";
push @usage, "  -t, --bedTools   Path to bed tools.\n";
push @usage, "  -a1, --UCSCAnnoFile text format file extracted from UCSC (UCSC, Tables, RefSeq, all fields from sekected table, txt output) .\n";
push @usage, "  -r1, --ISFile csv File with chr,IS, strand, seq count (sample file, the file generated by GENIS RescultsClusteredAnnotated.csv).\n\n";
push @usage, "This program takes as input the insertion site file in bed format and annotate them.\n\n";
push @usage, "#########################################################################\n\n";

my $help;
my $output_dir;
my $scriptDir;
my $libDir;
my $UCSCAnnoFile;
my $ISFile;

GetOptions
(
 'h|help|?'    => \$help,
 'o=s'    => \$output_dir,
 's=s'    => \$scriptDir,
 'l=s'    => \$libDir,
 't=s'    => \$bedTools,
 'a1=s'    => \$UCSCAnnoFile,
 'r1=s'    => \$ISFile,
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

if (defined $bedTools) {
        if (!-e $bedTools) {
                print "\nThe script directory $bedTools does not exist!\n\n";
                print @usage;
                exit;
        }
}else{
        print "Please provide path to bed tools !\n\n";
        print @usage;
        exit;
}


if (defined $UCSCAnnoFile) {
        if (!-e $UCSCAnnoFile){
                print "\nThe file $UCSCAnnoFile does not exist!\n\n";
                print @usage;
                exit;
        }
}
if (defined $ISFile) {
        if (!-e $ISFile){
                print "\nThe file $ISFile does not exist!\n\n";
                print @usage;
                exit;
        }
}
###############################################################################################
#print "\nIS annotation in progress...\n" ;
#print "\nResults will be stored in  $output_dir\n";

#modifying UCSC text table file to a bed file per requirement
system (" cut -f3,5,6 $UCSCAnnoFile > $output_dir/UCSCTableMod1.bed ");
system (" cut -f2,4,7,8,9,10,11,13 $UCSCAnnoFile > $output_dir/UCSCTableMod2.bed ");
system (" paste $output_dir/UCSCTableMod1.bed   $output_dir/UCSCTableMod2.bed  >  $output_dir/UCSCTableMod3.bed ");
system (" sed 's/ \+ /\t/g'  $output_dir/UCSCTableMod3.bed  > $output_dir/UCSCTableMod4a.bed ");
system (" sed  '1d'  $output_dir/UCSCTableMod4a.bed >  $output_dir/UCSCTableMod4.bed  ");

#modifying IS File to a bed file per requirement
system (" sed 's/^ *//g'  $ISFile | sed  's/ /\t/g'  > $output_dir/ISFileMod0_bed ");
system (" sed  's/ /\t/g' $ISFile| awk '{print (\$1,\$2, \$2+1, \$3, \$4) }' |  sed  's/ /\t/g' > $output_dir/ISFileMod1.bed ");
system (" $bedTools closest -a $output_dir/ISFileMod1.bed -b $output_dir/UCSCTableMod4.bed  -t first > $output_dir/anno1a.bed ");
system (" sed 's/ /\t/g' $output_dir/anno1a.bed > $output_dir/anno1.bed ");
#gene-length
system (" awk '{print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$8-\$7) }' $output_dir/anno1.bed > $output_dir/anno2.bed ");
system (" sed 's/ /\t/g' $output_dir/anno2.bed > $output_dir/anno3.bed ");
#distToTSS
system (" awk '{ if (\$2>=\$7 && \$2<=\$8 && \$10 ~ /+/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$2-\$7-1) ; if (\$2>=\$7 && \$2<=\$8 && \$10 ~ /-/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$8-\$2); if (\$2<\$7 || \$2>\$8) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17, \$10 )} '  $output_dir/anno3.bed > $output_dir/anno4.bed");
system (" sed 's/ /\t/g' $output_dir/anno4.bed > $output_dir/anno5.bed ");
#Upstream and Downstream distance
system ("awk '{ if ( \$2<\$7 && \$10  ~ /+/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$7-\$2+1); if (\$2>\$8 && \$10  ~ /-/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$2-\$8) ; if ( \$2>=\$7 && \$2<=\$8 && \$10 ~ /+/ || \$2>=\$7 && \$2<=\$8 && \$10 ~ /-/ ) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$10); if ( \$2>\$8 && \$10 ~ /+/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$10); if ( \$2<\$7 && \$10~ /-/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$10)} '  $output_dir/anno5.bed  > $output_dir/anno6.bed ");
system (" sed 's/ /\t/g' $output_dir/anno6.bed > $output_dir/anno7.bed ");

system ("awk '{ if ( \$2>\$8 && \$10   ~ /+/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$19 ,\$2-\$8); if ( \$2<\$7 && \$10 ~ /-/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$19 ,\$7-\$2+1); if ( \$2>=\$7 && \$2<=\$8 && \$10 ~ /+/ || \$2>=\$7 && \$2<=\$8 && \$10 ~ /-/  ) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$19 ,  \$10) ; if ( \$2<\$7 && \$10  ~ /+/) print (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$19,  \$10) ; if (\$2>\$8 && \$10  ~ /-/) print  (\$1,\$2, \$3, \$4 , \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13, \$14, \$15, \$16, \$17,\$18,\$19 ,  \$10)  } ' $output_dir/anno7.bed > $output_dir/anno8.bed ");
system (" sed 's/ /\t/g' $output_dir/anno8.bed > $output_dir/anno9a.bed ");

system(" sort -k1,1 -k2,2n  -k3,3 -k4,4n   $output_dir/ISFileMod0_bed >  $output_dir/ISFileMod00.bed");
system(" cut -f7- $output_dir/ISFileMod00.bed > $output_dir/ISFileMod000.bed");

system(" sed  's/\t//g' $output_dir/ISFileMod000.bed >  $output_dir/ISFileMod0000.bed");
system(" sort -k1,1 -k2,2n  -k3,3 -k4,4n  $output_dir/anno9a.bed > $output_dir/anno9b.bed");
system(" paste  $output_dir/anno9b.bed $output_dir/ISFileMod0000.bed > $output_dir/anno9.bed ");
#to assign exon and intron
system(" bash $scriptDir/exon.awk  $output_dir  anno9.bed ");
system ("  $bedTools intersect  -a $output_dir/exonsmyFile.bed  -b   $output_dir/anno9.bed  -wb  > $output_dir/exonsmyFile1a.bed  ");
system (" sort -k7,7 -k8,8n -k9,9n -k10,10 -k11,11 -k27,27  -u $output_dir/exonsmyFile1a.bed > $output_dir/exonsmyFile1.bed   ");
system(" bash $scriptDir/intron.awk  $output_dir  anno9.bed ");
system ("  $bedTools intersect  -a  $output_dir/intronsmyFile.bed  -b  $output_dir/anno9.bed -wb > $output_dir/intronsmyFile1a.bed  ");
system (" sort -k7,7 -k8,8n -k9,9n -k10,10 -k11,11 -k27,27  -u $output_dir/intronsmyFile1a.bed  > $output_dir/intronsmyFile1.bed   ");
system (" cut -f7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 $output_dir/exonsmyFile1.bed >  $output_dir/exonsmyFile11.bed ");
system (" cut -f1,2,3,4,5,6 $output_dir/exonsmyFile1.bed >  $output_dir/exonsmyFile11a.bed ");
system (" paste $output_dir/exonsmyFile11.bed $output_dir/exonsmyFile11a.bed  > $output_dir/exonsmyFile11b.bed  ");
system (" cut -f7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 $output_dir/intronsmyFile1.bed > $output_dir/intronsmyFile11.bed  ");
system (" cut -f1,2,3,4,5,6 $output_dir/intronsmyFile1.bed >  $output_dir/intronsmyFile11a.bed ");
system (" paste $output_dir/intronsmyFile11.bed $output_dir/intronsmyFile11a.bed  > $output_dir/intronsmyFile11b.bed  ");
system (" cut -f2  $output_dir/exonsmyFile11b.bed >  $output_dir/testt1 ");
system (" cut -f2  $output_dir/intronsmyFile11b.bed >  $output_dir/testt2 ");
system (" grep -v -Fwf $output_dir/testt1 $output_dir/testt2  >  $output_dir/testt3 ");
system (" awk 'NR==FNR{tgts[\$1]; next} \$2 in tgts' $output_dir/testt3  $output_dir/intronsmyFile11b.bed > $output_dir/intronsmyFile11b.ed.bed  ");
system (" mv $output_dir/intronsmyFile11b.ed.bed  $output_dir/intronsmyFile11b.bed ");
system (" cat $output_dir/exonsmyFile11b.bed  $output_dir/intronsmyFile11b.bed >  $output_dir/intronsExonsMyFile.bed ");
system (" echo notApplicable > $output_dir/temp.txt ");
system (" cat  $output_dir/intronsExonsMyFile.bed  $output_dir/temp.txt >  $output_dir/intronsExonsMyFile1.bed ");
system (" awk -F' ' 'NR==FNR{c[\$1\$2]++;next};c[\$1\$2] == 0' $output_dir/intronsExonsMyFile1.bed $output_dir/anno9.bed >  $output_dir/notInGeneMyFile.bed ");
system (" cat $output_dir/intronsExonsMyFile.bed  $output_dir/notInGeneMyFile.bed  > $output_dir/anno10.bed ");
system (" sort -k1,1 -k2,2  $output_dir/anno10.bed   > $output_dir/anno11.bed ");
system (" cut -f1,2,4,5,9,10,16,17,18,19,20,24  $output_dir/anno11.bed  >  $output_dir/anno12.bed  ");
system (" bash $scriptDir/formatAnnotation.sh $output_dir anno12.bed ISFileMod0_bed");
####################################################################################################



