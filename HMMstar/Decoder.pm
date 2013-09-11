package HMMstar::Decoder;
use strict;
use warnings 'FATAL' => 'all';

my $MIN_PROB  = 0;
my $MIN_SCORE = -999;

sub new {
	my ($class, $hmm, $dna) = @_;
	
	# re-map initial and terminal probabilities and transition matrix
	my %init;
	my %term;
	my %next;
	my %prev;
	foreach my $from (keys %{$hmm->{tmatrix}}) {
		foreach my $to (keys %{$hmm->{tmatrix}{$from}}) {
			if ($from eq 'INIT') {
				$init{$to} = $hmm->{tmatrix}{'INIT'}{$to};
			} elsif ($to eq 'TERM') {
				$term{$from} = $hmm->{tmatrix}{$from}{'TERM'};
			} else {
				$next{$from}{$to} = $hmm->{tmatrix}{$from}{$to};
				$prev{$to}{$from} = $hmm->{tmatrix}{$from}{$to};
			}
		}
	}
	
	# re-map states
	my %state;
	foreach my $s (@{$hmm->{states}}) {
		$state{$s->{name}} = $s;
	}
	
	my $self = bless {
		hmm     => $hmm,
		dna     => $dna,
		init    => \%init,
		term    => \%term,
		tnext   => \%next,
		tprev   => \%prev,
		state   => \%state,
		trellis => [],
		logs    => 0,
	};
	
	return $self;
}

sub decode_prob {
	my ($self) = @_;
	
	viterbi_prob($self);
	# forward
	# backward
	# posterior
}

sub decode_log {
	# viterbi
	# forward
	# backward
	# posterior
}

sub viterbi_prob {
	my ($self) = @_; 
	
	my $t = $self->{trellis};
	
	# initialization
	foreach my $s (keys %{$self->{state}}) {
		if (defined $self->{init}{$s}) {
			$t->[0]{$s}{vprob}  = $self->{init}{$s};
			$t->[0]{$s}{vtrace} = $s;
		} else {
			$t->[0]{$s}{vprob} = $MIN_PROB;
		}
	}
	
	# induction
	for (my $i = 1; $i <= $self->{dna}->length; $i++) {
		foreach my $s (keys %{$self->{state}}) {
			my ($prob, $trace) = viterbi_prob_calc($self, $i, $s);
			$t->[$i]{$s}{vprob} = $prob;
			$t->[$i]{$s}{vtrace} = $trace;
		}
	}
	
	# termination
	my $last = @$t;
	foreach my $s (keys %{$self->{state}}) {
		if (defined $self->{term}{$s}) {
			$t->[$last]{$s}{vprob} = $t->[$last-1]{$s}{vprob} * $self->{term}{$s};
		} else {
			$t->[$last]{$s}{vprob} = 0;
		}
	}
	
	# trace back
}	

sub viterbi_prob_calc {
	my ($self, $pos, $state) = @_;
	
	die "viterbi proc $pos $state\n";
	
	# emission probability
	my $emodel_name = $self->{state}{$state}{name};
	my $emodel = $self->{hmm}{emit}{$emodel_name};
	my $econtext = $self->{hmm}{context}{$emodel_name};
	
	my $nt = substr($self->{dna}{sequence}, $pos -1, 1);
	my $eprob;
	if ($econtext == 0) { 
		$eprob = $emodel->{""}{$nt};
	} elsif ($pos <= $econtext) {
		$eprob = 0.25; # if not enough context, use ambiguous value
	} else {
		my $ctx = substr($self->{dna}{sequence}, $pos - $econtext -1, $econtext);
		$eprob = $emodel->{$ctx}{$nt};
	}
	
	#print "pos = $pos, state = $state, nt = $nt, eprob = $eprob\n";
	
	# find max
	my $max_prob = 0;
	my $max_state = "";
	
	
	foreach my $from (keys %{$self->{tprev}}) {
		my $total_prob = $eprob * $self->{tprev}{$from}{$state} * 
			$self->{trellis}[$pos -1]{$from}{vprob};
		if ($total_prob > $max_prob) {
			$max_prob = $total_prob;
			$max_state = $from;
		}
	#	print "$from tprob = $total_prob\n";
	}
	#print "max: $max_state $max_prob\n";
	
	return $max_prob, $max_state;
}

sub debug {
	my ($self) = @_;
	
	my @t = @{$self->{trellis}};
	for (my $i = 0; $i < @t; $i++) {
		print "$i\t";
		if ($i == 0) {print "init\t"}
		elsif ($i == @t -1) {print "term\t"}
		else         {print substr($self->{dna}{sequence}, $i -1, 1), "\t"}
		foreach my $s (sort keys %{$t[$i]}) {
			printf "%.2g %s\t", $t[$i]{$s}{vprob}, $t[$i]{$s}{vtrace};
		}
		print "\n";
	}
}


1;
