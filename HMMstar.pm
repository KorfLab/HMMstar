###############################################################################
# HMMstar.pm
#
# ()_()
# (-.-)  Copyright 2007 Ian Korf. All rights reserved.
# (> <),
# 
#
###############################################################################
package HMMstar;

use strict;
use warnings 'FATAL' => 'all';

use HMMstar::DNA;
use HMMstar::Feature;
use HMMstar::Genome;
use HMMstar::Contig;
use HMMstar::Transcript;
use HMMstar::CodingTranscript;
use HMMstar::CDS;

use HMMstar::HMM;
use HMMstar::State;
use HMMstar::SequenceModel;
use HMMstar::Decoder;

###############################################################################

# The base class has no constructors, just general utility functions

#################
# DNA utilities #
#################

sub reverse_complement {
	my ($seq) = @_;
	$seq = reverse $seq;
	$seq =~ tr[ACGTRYMKWSBDHVacgtrymkwsdbhv]
	          [TGCAYRKMWSVHDBtgcayrkmwsvhdb];
	return $seq;
}

sub dna_word_table {
	my ($mer, $pseudo) = @_;
	
	$pseudo = 0 unless defined $pseudo;
	
	my %H;
	my $code = "";
	my @var;
	for (my $i = 0; $i < $mer; $i++) {
		my $tab = "\t" x $i;
		push @var, "\$c$i";
		$code .= "$tab foreach my \$c$i qw(A C G T) {\n";
	}
	$code .= "\t" x ($mer);
	my $var = join("", @var);
	$code .= "\$H{\"$var\"} = $pseudo\n";
	$code .= "}" x $mer;
		
	eval $code;
	return \%H;
}


#########################
# Translation utilities #
#########################

# Translation uses the standard genetic code. To use other genetic codes,
# use the edit_translation() function as needed. I should make this easier.

my %Translation = (
	'AAA' => 'K', 'AAC' => 'N', 'AAG' => 'K', 'AAT' => 'N', 
	'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T', 
	'AGA' => 'R', 'AGC' => 'S', 'AGG' => 'R', 'AGT' => 'S', 
	'ATA' => 'I', 'ATC' => 'I', 'ATG' => 'M', 'ATT' => 'I', 
	'CAA' => 'Q', 'CAC' => 'H', 'CAG' => 'Q', 'CAT' => 'H', 
	'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P', 
	'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 
	'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 
	'GAA' => 'E', 'GAC' => 'D', 'GAG' => 'E', 'GAT' => 'D', 
	'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A', 
	'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G', 
	'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V', 
	'TAA' => '*', 'TAC' => 'Y', 'TAG' => '*', 'TAT' => 'Y', 
	'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 
	'TGA' => '*', 'TGC' => 'C', 'TGG' => 'W', 'TGT' => 'C', 
	'TTA' => 'L', 'TTC' => 'F', 'TTG' => 'L', 'TTT' => 'F'
);

sub edit_translation {
	my ($codon, $aa) = @_;
	die "not a codon" unless $codon =~ /^[ACGT][ACGT][ACGT]$/;
	die "not an amino acid" unless $aa =~ /^[ACDEFGHIKLMNPQRSTVWY]$/;
	$Translation{$codon} = $aa;
}

sub translate {
	my ($seq, $start) = @_;
	
	my $trans = "";
	for (my $i = $start; $i < length($seq); $i+=3) {
		my $codon = uc(substr($seq, $i, 3));
		last if length($codon) < 3;
		if (not exists $Translation{$codon}) {$trans .= 'X'}
		else                                 {$trans .= $Translation{$codon}}
	}
	return $trans;
}

############################
# Error checking utilities #
############################

sub is_similar {
	my ($float1, $float2) = @_;
	if (abs($float1 - $float2) < 0.001) {return 1}
	return 0;
}

sub is_probability {
	my ($val) = @_;
	if ($val >= 0 and $val <= 1) {return 1}
	else                         {return 0}
}

sub is_integer {
	my ($val) = @_;
	if ($val == int($val)) {return 1}
	else                   {return 0}
}

1;


=head1 Name

HMMstar

=head1 Description

HMMstar is an HMM-based sequence analysis software package designed for
genome analysis. This module contains a variety of functions and classes
for building traditional and generalized HMMs.

=head1 Design Principles

HMMstar is designed to be (1) useful (2) educational (3) extensible.

=head1 Examples

=cut


