#!/usr/bin/perl
use strict;
use warnings;
use HMMstar;
use Getopt::Std;
our ($opt_h, $opt_v);
getopts('hv');

die "
usage: $0 [options] <HMMstar hmm file> <fasta file>
options:
  -h  help
  -v  verbose
" unless @ARGV == 2;

my ($HMM_FILE, $DNA_FILE) = @ARGV;

my $hmm = new HMMstar::HMM($HMM_FILE);
#my $dna = new HMMstar::DNA($DNA_FILE);
#my $decoder = new HMMstar::Decoder($hmm, $dna);

__END__
