package HMMstar::Transcript;
use strict;
use warnings 'FATAL' => 'all';

# Transcript object is a collection of exon feature objects. The Transcript object
# is created using the exon definition. Transcript object does NOT have to be coding. 

sub new {
	
	my ($class, $unsorted_exons) = @_;
	
	# exons must be sorted from this point
	my @exon = sort {$a->{start} <=> $b->{start}} @$unsorted_exons;
	
	my $self = bless {
		name => $exon[0]{group},
		strand => $exon[0]{strand},
		source => $exon[0]{source},
		exons  => \@exon, 	
		tss => $exon[0]{start},				# tss: transcription start site
		tts => $exon[0]{start},				# tts: transcription termination site
	};
 
 
 	my @tmp = reverse(@exon);
 	$self->{tss} = $tmp[0]{end} if ($self->{strand} eq '-');
 	$self->{tts} = $tmp[0]{end} if ($self->{strand} eq '+');
			
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

return $self;			
}


sub validate {
	my ($self) = @_;
	
	# check for introns that are way too short
	# do a lot of checks
}

sub exons {
	my ($self) = @_;
	return @{$self->{exons}};

}

sub name {
	my ($self) = @_;
	return $self->{name};
}

sub strand {
	my ($self) = @_;
	return $self->{strand};
}

sub start {
	my ($self) = @_;
	return $self->{tss};
}

sub end {
	my ($self) = @_;
	return $self->{tts};
}


1;
