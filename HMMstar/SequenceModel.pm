package HMMstar::SequenceModel;
use strict;
use warnings 'FATAL' => 'all';

sub new {
	my ($class, $name, $type, $order, $seqs, $pseudo) = @_;
	my $self = bless {
		name  => $name,
		type  => $type,
		order => $order,
	};
	
	if ($type eq 'MM') {
		$self->{model} = markov_model($order, $seqs, $pseudo);
	} elsif ($type eq 'PWM') {
		$self->{model} = pwm_model($order, $seqs, $pseudo);
	} else {
		die "$type unsupported at this time";
	}
	
	return $self;
}

sub write {
	my ($self, $file) = @_;

	open(OUT, '>', $file) or die;

	for my $prefix ( sort keys %{$self->{model}} ){
	    for my $nt ( sort keys %{$self->{model}->{$prefix}}  ) {
		printf OUT ("%s %s %.3f\n", $prefix, $nt, $self->{model}->{$prefix}{$nt});
	    }
	}
	close(OUT);
		

}

sub markov_model {
	my ($order, $seqs, $pseudo) = @_;
	$pseudo = 0 unless defined $pseudo;
	my @alph = qw(A C G T);
	my $table = HMMstar::dna_word_table($order, "{}");
	foreach my $prefix (keys %$table) {
		foreach my $nt (@alph) {
			$table->{$prefix}{$nt} = $pseudo;
		}
	}
	
	foreach my $seq (@$seqs) {
		for (my $i = 0; $i < length($seq) - $order -1; $i++) {
			my $prefix = substr($seq, $i, $order);
			my $nt = substr($seq, $i + $order + 1, 1);
			next unless $prefix =~ /^[ACGT]+$/;
			next unless $nt =~ /^[ACGT]$/;
			$table->{$prefix}{$nt}++;
		}
	}
	
	foreach my $prefix (keys %$table) {
		my $total = 0;
		my $zerocount = 0;
		foreach my $nt (@alph) {
			if ($table->{$prefix}{$nt} == 0) {
				$zerocount = 1;
			} else {
				$total += $table->{$prefix}{$nt};
			}
		}
		if ($zerocount) {
			die "some values in Markov model with zero counts, use pseudocounts";
		}
		
		foreach my $nt (@alph) {
			$table->{$prefix}{$nt} /= $total;
		}
	}
	
	
		
	return $table;
}


1;
