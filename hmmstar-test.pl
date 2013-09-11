#!/usr/bin/perl
use strict;
use warnings;
use HMMstar;

my $genome = new HMMstar::Genome("At.dna", "At.gff");


__END__
my $s1 = new HMMstar::State(
	'coding',
	0.5,
	0.5,
	{'coding' => 0.99, 'non-coding' => 0.01},
	{'' => [0.2, 0.3, 0.3, 0.2]},
);

my $s2 = new HMMstar::State(
	'non-coding',
	0.5,
	0.5,
	{'non-coding' => 0.99, 'coding' => 0.01},
	{'' => [0.3, 0.2, 0.2, 0.3]},
);

my $hmm = new HMMstar::HMM(
	'coding-2-0',
	'Ian Korf',
	'Simplest 2-state coding/non-coding',
	[$s1, $s2],
);

$hmm->write(\*STDOUT);

__END__
	my ($class, $name, $init, $term, $tran, $emit, $durn) = @_;
