# README #
(GENE-IS_1.0)

GENE-IS is a pipeline for the extraction of integration sites from next-generation sequencing data of clinical and preclinical gene therapy studies.
It is specifically designed in order to accept the sequencing reads originated from different protocols like LAM (linear amplification mediated) PCR and Targeted Sequencing (SureSelect/AGILENT) methods.


## How do I get set up? ###

### Installation

The easiest way to obtain and run gene-is is cloning the present repository
```
mkdir path_to_location
cd path_to_location
git clone https://github.com/G100DKFZ/gene-is.git
cd gene-is
```


## Testing

In order to test GENE-IS installation was successful;

* Type the following command on terminal for changing directory to scripts
```
cd /path_to_location/gene-is/scripts
# export the location of gene-is
export GENIS=/path_to_location/gene-is
# Run test suite by following command
./testGenis.sh
```
On the terminal will appear these options;
```
1) Targeted Sequencing Pair BWA 4) All
2) Targeted Sequencing Single 5) Clear
3) LAM-PCR 6) Quit
```

* To run tests for targeted sequencing paired end mode type at terminal 1 and press enter.
If installation was successful following message will appear on the terminal
"Targeted Sequencing Pair worked as expected!"



* To run tests for targeted sequencing single end mode type at terminal 2 and press enter.
If installation was successful following message will appear on the terminal
"Targeted Sequencing Single end worked as expected!"



* To run tests for LAM-PCR paired end mode type at terminal 3 and press enter.
If installation was successful following message will appear on the terminal
"LAM-PCR Pair worked as expected!"



## Dependencies

### Third-party tools

GENE-IS depends on several third party tools which are open source and are freely available.
All these tools are already provided within the GENE-IS package in tools/bin directory.

For user information names of tools and related links are provided here;
Tool 		Version 	URL

BWA 		0.7.4 		http://sourceforge.net/projects/bio-bwa/files/?source=navbar

Bedtools	2.17.0 	https://code.google.com/p/bedtools/downloads/detail?name=BEDTools.v2.17.0.tar.gz&can=2&q=

Samtools 	0.1.19 	http://samtools.sourceforge.net/

BLAT 		v.35 		http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip

Skewer		0.1.117 	http://sourceforge.net/projects/skewer/files/Binaries/



### Perl modules

The required Perl libraries are pre-packaged within the tool ("lib" dir in GENE-IS).


## Configuration File
GENE-IS has specific configuration files for each mode of analysis; LAM-PCR, TES paired and TES single end configuration files.
Only the relevant configuration file should be modified for particular analysis.
For testing GENE-IS installation user does not need to change any paramter in the configurtaion file.
The templates are in the gene-is path; i.e.
```
$GENIS/configFile_targetedSequencing_pairedEnd.txt
```
* Analysis

> In order to test other datasets used in the Manuscript for GENE-IS benchmarking, please see "README" file in "/path_to_location/gene-is/testFiles" directory


### Synthetic and test guidelines

* Synthetic datasets
* Other guidelines

### Contacts ###

* Contact: raffaele.fronza@nct-heidelberg.de
* Contact: saira.afzal@nct-heidelberg.de
