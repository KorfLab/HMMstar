package HMMstar::CDS;
use strict;
use warnings 'FATAL' => 'all';

sub best_orf {
	die "best_orf not implemented yet";
}

sub new {
	my ($class, $unsorted_exons) = @_;
	
	# exons must be sorted from this point
	my @exon = sort {$a->{start} <=> $b->{start}} @$unsorted_exons;
	
	my $self = bless {
		name => $exon[0]{group},
		strand => $exon[0]{strand},
		source => $exon[0]{source},
		exons  => \@exon,
	};
	
	
	# check for catastrophic errors, must abort if found
	my $mixed_strand = 0;
	my $overlap = 0;
	for (my $i = 1; $i < @exon; $i++) {
		if ($exon[0]{strand} ne $exon[$i]{strand}) {$mixed_strand = 1}
		if (HMMstar::Feature::overlap($exon[$i-1], $exon[$i])) {$overlap = 1}
	}
	if ($mixed_strand or $overlap) {
		$self->{error} = 1;
		if ($mixed_strand) {warn "mixed strand in $exon[0]{group}"}
		if ($overlap) {warn "overlapping exons in $exon[0]{group}"}
		return $self;
	}
		
	# create introns
	my @intron;
	for (my $i = 1; $i < @exon; $i++) {
 
		my $ib = $exon[$i-1]{end} + 1;	# $ib is the intron begin index
		my $ie = $exon[$i]{start} -1;	# $ie is the intron end index		
		push @intron, new HMMstar::Feature($self->{source}, 'Intron',
			$ib, $ie, $self->{strand}, 0, $self->{name});
	}
	$self->{introns} = \@intron;
	
	# transcribe & translate
	my $tx = "";

	if ($self->{strand} eq '+') {
		foreach my $exon (@{$self->{exons}}) {$tx = $tx . $exon->sequence}
	}

	if ($self->{strand} eq '-') {
		foreach my $exon (@{$self->{exons}}) {$tx = $exon->sequence . $tx}
	}
	
	$self->{transcript} = $tx;
	$self->{protein} = HMMstar::translate($tx, 0); # will validate later
	
	
	$self->validate;
	return $self;
	
}

sub protein {
	my ($self) = @_;
	return $self->{protein};
}

sub transcript {
	my ($self) = @_;
	return $self->{transcript};
}

sub name {
	my ($self) = @_;
	return $self->{name};
}

sub strand {
	my ($self) = @_;
	return $self->{strand};
}

sub validate {
	my ($self) = @_;
	
	# check for introns that are way too short
	# do a lot of checks
	
	# This will be done by a driver script that uses various objects of HMMstar to do extensive 
	# validation of a given potential CDS
}

# gene features

sub introns {
	return @{shift->{introns}};
}

sub exons {
	return @{shift->{exons}};
}

sub start_site {
	my ($self) = @_;
	
	# start site defined as the 'A' of 'ATG'
	my $pos;
	if ($self->{strand} eq '+') {
		$pos = $self->{exons}[0]{start};
	} else {
		$pos = $self->{exons}[@{$self->{exons}}-1]{end};
	}
	
	my $site = new HMMstar::Feature(
		$self->{source},
		'Start',
		$pos,
		$pos,
		$self->{strand},
		0,
		$self->{group});
	return $site;
	
}

sub stop_site {
	my ($self) = @_;
	
	# stop site defined as the first 'T' in {TAA, TGA, TAG}
	my $pos;
	if ($self->{strand} eq '+') {
		$pos = $self->{exons}[@{$self->{exons}}-1]{end} -2;
	} else {
		$pos = $self->{exons}[0]{start} +2;
	}
	
	my $site = new HMMstar::Feature(
		$self->{source},
		'Stop',
		$pos,
		$pos,
		$self->{strand},
		0,
		$self->{group});
	return $site;
}

sub acceptor_sites {
	my ($self) = @_;
	
	# acceptor defined as the 'G' in the '...AG' consensus
	my @f;
	foreach my $intron (@{$self->{introns}}) {
		my $pos;
		if ($self->{strand} eq '+') {
			$pos = $intron->{end};
		} else {
			$pos = $intron->{start};
		}
		push @f, new HMMstar::Feature(
			$self->{source},
			'Acceptor',
			$pos,
			$pos,
			$self->{strand},
			0,
			$self->{group});
	}
	
	@f = reverse @f if $self->{strand} eq '-';
	return @f;
}

sub donor_sites {
	my ($self) = @_;
	
	# donor site defined as the 'G' in the 'GT...' consensus
	my @f;
	foreach my $intron (@{$self->{introns}}) {
		my $pos;
		if ($self->{strand} eq '+') {
			$pos = $intron->{start};
		} else {
			$pos = $intron->{end};
		}
		push @f, new HMMstar::Feature(
			$self->{source},
			'Donor',
			$pos,
			$pos,
			$self->{strand},
			0,
			$self->{group});
	}
	
	@f = reverse @f if $self->{strand} eq '-';
	return @f;
}


1;
