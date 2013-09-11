package HMMstar::DNA;
use strict;
use warnings 'FATAL' => 'all';

sub new {
	my ($class, $id, $seq) = @_;
	my $self = bless {
		definition => $id,
		sequence => $seq,
		'length' => CORE::length($seq),
	};
	$self->validate;
	return $self;
}

sub read {
	my ($FH) = @_;
	my $head = <$FH>;
	my ($def) = $head =~ /^>(.+)/;
	die "improper FASTA format" unless defined $def;
	my $seq = "";
	while (<$FH>) {
		chomp;
		if (/^[ACTGRYMKWSBDHVNactgrymkwsbdhvn]+$/) {$seq .= $_}
		else {die "HMMstar::DNA::read error $_"}
	}	
	return new HMMstar::DNA($def, $seq);
}

sub validate {
	my ($self) = @_;
	
	if ($self->{sequence} !~ /^[ACTGRYMKWSBDHVNactgrymkwsbdhvn]+$/) {
		die "sequence has non-DNA symbols";
	}
	if ($self->{definition} !~ /^\S+/) {
		die "sequence has no definition";
	}
}

sub definition {
	my ($self, $val) = @_;
	return $self->{definition} unless defined $val;
	$self->{definition} = $val;
	$self->validate;
}

sub sequence {
	my ($self, $val) = @_;
	return $self->{sequence} unless defined $val;
	$self->{sequence} = $val;
	$self->validate;
}

sub length {
	my ($self) = @_;
	return $self->{'length'}
}

sub reverse_complement {
	my ($self) = @_;
	my $rc = HMMstar::reverse_complement($self->{sequence});
	return new HMMstar::DNA($self->{definition}, $rc);
}

sub translate {
	my ($self, $offset) = @_;
	return HMMstar::translate($self->{sequence}, $offset);
}

1;
