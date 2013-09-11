package HMMstar::HMM;
use strict;
use warnings 'FATAL' => 'all';

sub read {
	my ($FH) = @_;
	
	my $self = bless {};
	
	while (<$FH>) {
		next if /^#/;
		next unless /\S/;
		last;
	}
	
	my $head = $_;
	if ($head !~ /^(HMMstar-HMM|^denada-HMM)/) {
		die "unknown format\n";
	}
	
	# read states
	my ($type, $name, $states) = split;
	$self->{type} = $type;
	$self->{name} = $name;
	$self->{states} = [];
	for (my $i = 0; $i < $states; $i++) {
		push @{$self->{states}}, HMMstar::State::read($FH);
	}
	
	# construct transition matrix
	my %tm;
	foreach my $state (@{$self->{states}}) {
		my $s1 = $state->{name};
		my $total = 0;
		foreach my $s2 (keys %{$state->{trans}}) {
			$tm{$s1}{$s2} = $state->{trans}{$s2};
			$total += $state->{trans}{$s2};
		}
		if (not HMMstar::is_similar(1, $total)) {
			HMMstar::browse($self);
			die "transitions ($total) do not sum close enough to 1.000";
		}
	}
	$self->{tmatrix} = \%tm;
		
	return $self;
}

sub write {
	my ($self, $FH) = @_;
	die "writing not yet supported";
}


1;
