package HMMstar::State;
use strict;
use warnings 'FATAL' => 'all';

sub new {
	my ($class, $name, $init, $term, $links, $order, $duration) = @_;
	die "constructor not finished";
}

sub read {
	my ($FH) = @_;
	
	while (<$FH>) {last if /^\S/}
	my $head = $_;
	if ($head !~ /^HMMstar-State|^denada-State/) {
		die "unsupported format";
	}
	my ($type, $name, $init, $term, $links, $order, $limit) = split;
	
	my $self = bless {
		name => $name,
		init => $init,
		term => $term,
		links => $links,
		order => $order,
		limit => $limit,
	};
	
	# read trans
	my $read = 1; # lines
	my %trans;
	while (<$FH>) {
		next unless /\S/;
		my ($state, $prob) = split;
		$trans{$state} = $prob;
		last if $read++ == $links;
	}
	$self->{trans} = \%trans;
	
	# read model
	my @model;
	$read = 4 ** ($order +1); # values
	while (<$FH>) {
		next unless /\S/;
		my @val = split;
		push @model, @val;
		last if @model == $read;
	}
	$self->{model} = HMMstar::dna_word_table($order +1, 0);
	foreach my $word (sort keys %{$self->{model}}) {
		$self->{model}{$word} = shift @model;
	}

	# read range
	my @range;
	if ($limit) {
		while (<$FH>) {
			next unless /\S/;
			my @val = split;
			push @range, @val;
			last if @range == $limit
		}
	}
	$self->{range} = @range;
		
	return $self;
	
}

sub write {
}


1;

