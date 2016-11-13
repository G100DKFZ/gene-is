# README #
(EXTERNAL VERSION GENE-IS_1.0)

GENE-IS is a pipeline for the extraction of IS from NGS clinical and preclinical gene therapy studies. 
It is specifically designed in order to accept reads originated from different protocols like LAM-PCR (linear amplification mediated) and Targeted Sequencing (SureSelect/AGILENT) methods.

#######################################################################

### How do I get set up? ###

* Installation

Download GENE-IS source code by pasting this command on Linux terminal;

hg clone https://TesterGIS@bitbucket.org/dkfzto/gene-is1.0

Password; TesterGIS

#######################################################################

* Dependencies

Third-party tools:

GENE-IS depends on several third party tools which are open source and are freely available. 
User needs to download and install them before running GENE-IS. Go to the links of individual tools, download and install these tools by following instructions in the related tool manual/web page.

Tool 		Version 	URL

BWA 		0.7.4 		http://sourceforge.net/projects/bio-bwa/files/?source=navbar

Bedtools	2.17.0 		https://code.google.com/p/bedtools/downloads/detail?name=BEDTools.v2.17.0.tar.gz&can=2&q=

Samtools 	0.1.19 		http://samtools.sourceforge.net/

BLAT 		v.35 		http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip

Skewer		0.1.117 	http://sourceforge.net/projects/skewer/files/Binaries/

Perl modules:

The required Perl libraries are pre-packaged within the tool ("lib" dir in GENE-IS).

Only two module, need to be installed by the user are; Bio::SeqIO and Bio::DB::Sam 

These are usually available as a package on the Linux distributions.

However, can also be obtained from these links;

Bio::DB::Sam: http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Sam.pm

Note: Please install SAMtools library before installing Bio::DB::Sam module.

Bio::SeqIO: http://search.cpan.org/~cjfields/BioPerl-1.6.924/Bio/SeqIO.pm
#######################################################################

* Configuration File
GENE-IS has specific configuration files for each mode of analysis; LAM-PCR, TES paired and TES single end configuration files. 

Only the relevant configuration file should be modified for particular analysis. 


To test LAM-PCR mode of analysis, open “configFile_LAM-PCR_pairedEnd.txt”

User is required to provide correct paths for variables in the configuration file

Provide path to the raw data files (these files are in directory "/path_to_location/gene-is1.0/test/LAM-PCR/")

Also provide path to the reference genome (human genome hg38) index files for BWA and BLAT aligner (Please see manual section 3.2 how to create reference index files)

In addition, provide path to above mentioned third-party tools. 

* How to run tests for LAM-PCR mode:

To test LAMPCR mode of GENE-IS we will use a dataset of 600 reads. Follow these instructions to complete testing process.

• Create a working/result directory for testing LAM by following this command.

mkdir /path_to_location/gene-is1.0/test/LAM-PCR/results

• Type following command on terminal for changing directory to scripts

cd /path_to_location/gene-is1.0/scripts

• Type following commands on terminal

export GENIS=/path_to_location/gene-is1.0

• Run test suite by following command

./testGenis.sh

On the terminal will appear these options;

1) Targeted Sequencing Pair BWA 4) All 

2) Targeted Sequencing Single 5) Clear 

3) LAM-PCR 6) Quit

• Type at terminal 3 and press enter. Analysis will start and when it is finished type 6 at terminal to exit.

• Go to results directory by following command;

cd /path_to_location/gene-is1.0/test/LAM-PCR/results/

Open file with suffix “.GeneralStatistics.txt”, if you see these following lines in your file, it means the pipeline is correctly installed and LAM-PCR mode is working correctly.
##################################################
GENERAL STATISTICS
##################################################
Number of raw read pairs
600

Number of filtered and trimmed read pairs
500

Number of unclustered Integartion Sites (sequence_count)
461

Number of clustered Integartion Sites
461

##################################################



To test targeted sequencing mode (TES) of analysis, open “configFile_targetedSequencing_pairedEnd.txt” 

User is required to provide correct paths for variables in the configuration file

Provide path to the raw data files (these files are in directory "/path_to_location/gene-is1.0/test/targetedSequencing/")

Also provide path to the reference genome + vectorSeq (hg38+vectorSeq) index files for BWA and BLAT aligner (Please see manual section 3.2 how to create reference index files)

(The "vectorSeq.fa" is the vector file that is already present in directory "/path_to_location/gene-is1.0/test/targetedSequencing/")
In addition, provide path to above mentioned third-party tools.

* How to run tests for targeted sequencing mode:

To test TES paired end mode of GENE-IS we will use a one sample dataset of 500 raw reads. Follow these instructions to complete testing process.

• Create a working/result directory for testing TES by following this command.

mkdir /path_to_location/gene-is1.0/test/targetedSequencing/results/pairedEnd/

• Type following command on terminal for changing directory to scripts

cd /path_to_location/gene-is1.0/scripts

• Type following commands on terminal

export GENIS=/path_to_location/gene-is1.0

On the terminal will appear these options;

1) Targeted Sequencing Pair BWA 4) All

2) Targeted Sequencing Single 5) Clear

3) LAM-PCR 6) Quit

• Type at terminal 1 and press enter. Analysis will start and when it is finished type 6 at terminal to exit.

• Go to results directory by following command

cd /path_to_location/gene-is1.0/test/targetedSequencing/results/pairedEnd

Open file with suffix “.GeneralStatistics.txt”, if you see following initial lines in file, it means the pipeline is correctly installed and targeted sequencing mode is working succssfully.
#####################################################
GENERAL STATISTICS
#####################################################
Number of raw read pairs
500

Number of filtered and trimmed read pairs
500

Number of correctly aligned vector-vector read pairs
0

Number of unclustered Integration Sites (sequence_count)
490

Number of clustered Integration Sites
490
#####################################################



### In order to test other datasets used in Manuscript for GENE-IS benchmarking, please see "README" file in "/path_to_location/gene-is1.0/testFiles" directory ###


### Synthetic and test guidelines ###

* Synthetic datasets
* Other guidelines

### Contacts ###

* Contact: raffaele.fronza@nct-heidelberg.de
* Contact: saira.afzal@nct-heidelberg.de
