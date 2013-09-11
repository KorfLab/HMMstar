package HMMstar::Feature;
use strict;
use warnings 'FATAL' => 'all';

sub overlap {
	my ($f1, $f2) = @_;
	my ($x1, $x2) = ($f1->{start}, $f1->{end});
	my ($y1, $y2) = ($f2->{start}, $f2->{end});
	return 1 if
		($x1 >= $y1 and $x1 <= $y2) or
		($x2 >= $y1 and $x2 <= $y2);
	return 0;
}

my %Type = (
	'CDS'      			=> 1,
	'Intron'   			=> 1,
	'Exon'     			=> 1,
	'Donor'    			=> 1,	# 'G' in 'GT'
	'Acceptor' 			=> 1,	# 'G' in 'AG'
	'Start'    			=> 1,	# 'A' in 'ATG'
	'Stop'     			=> 1,	# 'T' in {'TAA,'TGA','TAG'}
	'FivePrimeExon' 	=> 1,
	'FivePrimeIntron'	=> 1,
	'ThreePrimeExon' 	=> 1,
	'ThreePrimeIntron'	=> 1,
);

sub new {
	my ($class, $source, $type, $start, $end, $strand, $score, $group) = @_;
		
	my $self = bless {
		source => $source,
		type   => $type,
		start  => $start,
		end    => $end,
		strand => $strand,
		score  => $score,
		group  => $group,
	};
		
	$self->validate;
	return $self;
}

sub validate {
	my ($self) = @_;
	
	my $max_coordinate = $self->{source}->length();	

	# check source
	die "source not HMMstar::DNA " unless ref $self->{source} eq 'HMMstar::DNA';
	
	# check type
	die "unknown type ($self->{type})" unless defined $Type{$self->{type}};
	
	# check start
	die "start not an integer" unless HMMstar::is_integer($self->{start});
	die "non-positive start coordinate ($self->{start})" if $self->{start} < 1;
	die "start out of range" if $self->{start} > $max_coordinate;

	# check end
	die "end not an integer" unless HMMstar::is_integer($self->{end});
	die "non-positive end coordinate ($self->{end})" if $self->{end} < 1;
	die "end out of range" if $self->{end} > $max_coordinate;

	# check strand
	die "unknown strand ($self->{strand})" unless $self->{strand} =~ /^[+-]$/;
	die "$self->{group} $self->{type}: end < start: $self->{end} < $self->{start}" if $self->{end} < $self->{start};
	
	# score is not checked now

	# group is optional and need not be checked

}

sub source {
	my ($self, $val) = @_;
	return $self->{source} unless defined $val;
	$self->{source} = $val;
	$self->validate;
}

sub type {
	my ($self, $val) = @_;
	return $self->{type} unless defined $val;
	$self->{type} = $val;
	$self->validate;
}

sub start {
	my ($self, $val) = @_;
	return $self->{start} unless defined $val;
	$self->{start} = $val;
	$self->validate;
}

sub end {
	my ($self, $val) = @_;
	return $self->{end} unless defined $val;
	$self->{end} = $val;
	$self->validate;
}

sub strand {
	my ($self, $val) = @_;
	return $self->{strand} unless defined $val;
	$self->{strand} = $val;
	$self->validate;
}

sub score {
	my ($self, $val) = @_;
	return $self->{score} unless defined $val;
	$self->{score} = $val;
	$self->validate;
}

sub group {
	my ($self, $val) = @_;
	return $self->{group} unless defined $val;
	$self->{group} = $val;
	$self->validate;
}

sub length {
	my ($self) = @_;
	return $self->{end} - $self->{start} + 1;
}

sub sequence {
	my ($self, $up, $down) = @_;
	$up   = 0 if not defined $up;
	$down = 0 if not defined $down;
	($up, $down) = ($down, $up) if $self->{strand} eq '-';
	my $subseq = substr(
		$self->{source}{sequence},
		$self->{start} -1 - $up,
		$self->{end} - $self->{start} +1 + $up + $down
	);
	if ($self->{strand} eq '-') {
		$subseq = HMMstar::reverse_complement($subseq);
	}
	return $subseq;
}


1;
