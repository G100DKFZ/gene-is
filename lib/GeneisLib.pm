#!/usr/bin/perl -w

package GeneisLib;

use diagnostics;
use strict;
use warnings;

$| = 1; #disable buffer

use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
our @ISA         = qw (Exporter);
our @EXPORT      = qw ($DEBUG_MODE new read has_value get_value average stdev Fastq2Fasta ChopReads ReformatOutput CalDustScore BlatOutCandidate BlastnCleanup PulloutVirus PulloutTop1Virus PulloutVirusList GetIdir GetChrPrefix GetInsertSize GetRefSubSeq FilterRefSubSeq WriteRefSubSeq ExtractMappedReads);



sub new {
  my $self = {};
  bless($self);
  return $self;
}


sub read {
  my $self = shift;
  my $config_file = shift;
  my %config_values;

  open CFG, $config_file or die "Error: Unable to open $config_file\n";
  while (<CFG>){
      chomp;
      /^\s*([^=\s]+)\s*=\s*(.*)$/;

      my $key = $1;
      my $value = $2;

      next if not defined $key;
      next if not defined $value;

      $config_values{$key} = $value;
  }
  close CFG;

  foreach my $key (keys %config_values){
      while ($config_values{$key} =~ /\$\(([^)]+)\)/){
	  my $other_key = $1;

	  if (not defined $config_values{$other_key}){
	      die "Error: no value for $other_key in config file $config_file\n";
	  }

	  $config_values{$key} =~ s/\$\($other_key\)/$config_values{$other_key}/;
      }
  }

  $self->{"config_values"} = \%config_values;
  $self->{"config_file"} = $config_file;
}

sub has_value {
  my $self = shift;
  my $key = shift;

  my $config_values = $self->{"config_values"};
  my $config_file = $self->{"config_file"};

  defined $config_values and defined $config_file or die "Error: config not read\n";

  return defined $config_values->{$key};
}

sub get_value {
  my $self = shift;
  my $key = shift;

  my $config_values = $self->{"config_values"};
  my $config_file = $self->{"config_file"};

  defined $config_values and defined $config_file or die "Error: config not read\n";

  return $config_values->{$key};
}

1;
