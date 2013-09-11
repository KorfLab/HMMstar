package HMMstar::Contig;
use strict;
use warnings 'FATAL' => 'all';

sub new {
	my ($class, $dna, $features) = @_;
	my $self = bless {
		dna => $dna,
		features => $features,
	};
	return $self;
}

sub dna {return shift->{dna}}
sub features {return @{shift->{features}}}

sub intergenic_sequences {
	my ($self) = @_;
	
	my $seq = uc $self->{dna}{sequence};
	
	# sort by group
	my %group;
	foreach my $feature (@{$self->{features}}) {
		push @{$group{$feature->{group}}}, $feature;
	}
	
	foreach my $gid (keys %group) {
		my ($min, $max) = (1e20, 0);
		foreach my $feature (@{$group{$gid}}) {
			$min = $feature->{start} if $feature->{start} < $min;
			$min = $feature->{end}   if $feature->{end}   < $min;
			$max = $feature->{start} if $feature->{start} > $max;
			$max = $feature->{end}   if $feature->{end}   > $max;
		}
		substr($seq, $min -1, $max - $min +1) = lc substr($seq, $min -1, $max - $min +1);
	}
	
	my @seq = split(/[a-z]+/, $seq);
	return @seq;
}

sub transcript_sequences {
	my ($self) = @_;
	
	# print "Inside transcript_sequences \n";
	
	# sort features by group name, keep only Coding_transcript features
	my %transcript;
	foreach my $feature (@{$self->{features}}) { 
		next unless ($feature->{type} eq 'Exon');
		push @{$transcript{$feature->{group}}}, $feature;
	}
		
	my @transcript;
	foreach my $group (keys %transcript) {
		push @transcript, new HMMstar::Transcript($transcript{$group});
	}
	
	return @transcript;
}

sub get_cds_info {
	my ($self) = @_;
	
	# get cds start stop coordinates 
	
	# sort features by group name, keep only CDS features
	my @cds;
	
	foreach my $feature (@{$self->{features}}) {
		next unless ($feature->{type} eq 'CDS');	
		push @cds, join(":",$feature->{group},$feature->{start},$feature->{end}) if ($feature->{strand} eq '+');
		push @cds, join(":",$feature->{group},$feature->{end},$feature->{start}) if ($feature->{strand} eq '-');
	}
			
	return @cds
}

sub coding_sequences {
	my ($self) = @_;
	
	# sort features by group name, keep only CDS features
	my %cds;
	foreach my $feature (@{$self->{features}}) {	
		next unless ($feature->{type} eq 'CDS');
		push @{$cds{$feature->{group}}}, $feature;
	}
		
	my @cds;
	foreach my $group (keys %cds) {
		push @cds, new HMMstar::CDS($cds{$group});
	}
	
	return @cds;
}

1;
