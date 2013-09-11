package HMMstar::CodingTranscript;
use strict;
use warnings 'FATAL' => 'all';

# The areas of the Transcript object that do not lie within the CDS are the 5' & 3' UTRs.

sub new {
	my ($class,$transcript,$cds2start,$cds2stop) = @_;	
	my @intron;
	
	my @AnnotatedExons;	# These exons have either Exon, ThreePrimeExon or FivePrimeExon types
	my @InternalExons;  # Include only exons in the CDS region
	my @UTR5; 			# Includes both exons and intron objects in 5'UTR 
	my @UTR3; 			# Includes both exons and intron objects in 3'UTR
	my $cdsObject;
	
	
	my $self = bless {
		transcript => $transcript,
		strand => $transcript->strand,
		cds2start => $cds2start,
		cds2stop => $cds2stop,
		introns => \@intron,
		tss => $transcript->start,
		tts => $transcript->end,
		cds => \$cdsObject,
		full => 0,
		utr3 => \@UTR3,
		utr5 => \@UTR5,
	};

	# Dump a lot of info to standard output
	print "\nINSIDE CodingTranscript CONSTRUCTOR \n";
	print "CDS start/end: $self->{cds2start},$self->{cds2stop} \n";
	print "STRAND = $self->{strand} \n";
	print "TSS: $self->{tss} \n";
	print "TTS: $self->{tts} \n";

	my @exons = $self->{transcript}->exons();
	my $i = -1;
	while (@exons) {
		++$i;
		my $exon = shift(@exons);
				
		print "Exon start/end: $exon->{start}, $exon->{end} \n";
		
		# Figure out if you are dealing with 5'UTR, CDS or 3'UTR feature exons. 
		# Treat each strand differently. If the CDS boundaries fall in the middle of the exon, 
		# break up that exon into two feature objects: one Exon feature object in the CDS region, 
		# one feature object either as FivePrimeExon or a ThreePrimeExon feature object.
 		   
		if ($self->{strand} eq '+') {
			if (($exon->{start} >= $self->{cds2start}) && ($exon->{end} <= $self->{cds2stop})){
				# Internal Exon
				push @AnnotatedExons, $exon; 
			}elsif (($exon->{start} < $self->{cds2start}) && ($exon->{end} < $self->{cds2start})) {
				# Exon lies entirely in the 5'UTR
				$exon->type('FivePrimeExon');
				push @AnnotatedExons, $exon;
			}elsif (($exon->{start} > $self->{cds2stop}) && ($exon->{end} > $self->{cds2stop})){
				# Exon lies entirely in the 3'UTR
				$exon->type('ThreePrimeExon');
				push @AnnotatedExons, $exon;
			}elsif (($exon->{start} < $self->{cds2start}) && ($exon->{end} > $self->{cds2start})){
				push @AnnotatedExons, new HMMstar::Feature($exon->{source}, 'FivePrimeExon',
				$exon->{start}, $self->{cds2start} - 1, '+', 0, $exon->{group});
				push @AnnotatedExons, new HMMstar::Feature($exon->{source}, 'Exon',
				$self->{cds2start}, $exon->{end}, '+', 0, $exon->{group});
			}elsif (($exon->{start} < $self->{cds2stop}) && ($exon->{end} > $self->{cds2stop})){
				push @AnnotatedExons, new HMMstar::Feature($exon->{source}, 'Exon',
				$exon->{start}, $self->{cds2stop}, '+', 0, $exon->{group});
				push @AnnotatedExons, new HMMstar::Feature($exon->{source}, 'ThreePrimeExon',
				$self->{cds2stop} + 1, $exon->{end}, '+', 0, $exon->{group});			
			}
		}
			
		if ($self->{strand} eq '-') {
			if (($exon->{start} >= $self->{cds2stop}) && ($exon->{end} <= $self->{cds2start})) {
				# Internal Exon  
				push @AnnotatedExons, $exon;
			}elsif (($exon->{start} > $self->{cds2stop}) && ($exon->{end} < $self->{cds2start})) {
				# Internal Exon
				push @AnnotatedExons, $exon;
			}elsif (($exon->{end} > $self->{cds2start}) && ($exon->{start} > $self->{cds2start})) {
				#print "Negative Case 1\n";
				$exon->type('FivePrimeExon');
				push @AnnotatedExons, $exon;			
			}elsif (($exon->{start} < $self->{cds2stop}) && ($exon->{end} < $self->{cds2stop})) {
				#print "Negative Case 2\n";
				$exon->type('ThreePrimeExon');
				push @AnnotatedExons, $exon;						
			} elsif (($exon->{start} < $self->{cds2stop}) && ($exon->{end} > $self->{cds2stop})) {
				#print "Negative Case 3\n";
				push @AnnotatedExons, new HMMstar::Feature($exon->{source}, 'ThreePrimeExon',
				$exon->{start}, $self->{cds2stop} + 1, '-', 0, $exon->{group});
				push @AnnotatedExons, new HMMstar::Feature($exon->{source}, 'Exon',
				$self->{cds2stop},$exon->{end}, '-', 0, $exon->{group});
			}elsif (($exon->{start} < $self->{cds2start}) && ($exon->{end} > $self->{cds2start})) {
				#print "Negative Case 4\n";
				push @AnnotatedExons, new HMMstar::Feature($exon->{source}, 'Exon',
				$exon->{start},$self->{cds2start}, '-', 0, $exon->{group});
				push @AnnotatedExons, new HMMstar::Feature($exon->{source}, 'FivePrimeExon',
				$self->{cds2start} + 1, $exon->{end},'-', 0, $exon->{group});
			}	
		}
	}
	
	
	print "\n\n";
	
	
	while (@AnnotatedExons) {
		my $exon = shift(@AnnotatedExons);
		print "$exon->{type} $exon->{start} $exon->{end} $exon->{strand} $exon->{group} \n";
		push @InternalExons, $exon if ($exon->{type} eq 'Exon');
		push @UTR5, $exon if ($exon->{type} eq 'FivePrimeExon');
		push @UTR3, $exon if ($exon->{type} eq 'ThreePrimeExon');
	}

	# if the coding transcript has both 5'UTR & 3'UTR then it is considered full
	$self->{full} = 1 if (@UTR3 && @UTR5);

	# Create a CDS object based on internal exons
	$cdsObject = new HMMstar::CDS(\@InternalExons);

	@exons = $self->{transcript}->exons();
	# Create introns in the 5UTR or 3UTR if there are any. Introns internal to the CDS will be 
	# created by the CDS object 
	for (my $i = 1; $i < @exons; $i++) {
		my $ib = $exons[$i-1]{end} + 1;				# Intron begin coordinate
		my $ie = $exons[$i]{start} -1;				# Intron end coordinate
		
		if ($self->{strand} eq '+') {
			if ($ib < $self->{cds2start}) {
				push @UTR5, new HMMstar::Feature($exons[$i]->source, 'FivePrimeIntron',
				$ib, $ie, $exons[$i]->strand, 0, $exons[$i]->group);
			}
			if ($ib > $self->{cds2stop}) {
				push @UTR5, new HMMstar::Feature($exons[$i]->source, 'ThreePrimeIntron',
				$ib, $ie, $exons[$i]->strand, 0, $exons[$i]->group);
			}						
		}
		
		if ($self->{strand} eq '-') {
			if ($ib > $self->{cds2start}) {
				push @UTR5, new HMMstar::Feature($exons[$i]->source, 'FivePrimeIntron',
				$ib, $ie, $exons[$i]->strand, 0, $exons[$i]->group);
			}
			if ($ib < $self->{cds2stop}) {
				push @UTR5, new HMMstar::Feature($exons[$i]->source, 'ThreePrimeIntron',
				$ib, $ie, $exons[$i]->strand, 0, $exons[$i]->group);
			}						
		}		
	}
	
	
	$self->validate;
	return $self;
	
}


# gene features

sub tss {
	my ($self) = @_;
	return $self->{tss};
}

sub tts {
	my ($self) = @_;
	return $self->{tts};
}


sub strand {
	my ($self) = @_;
	return $self->{strand};
}


sub get5UTR {
	my ($self) = @_;
	return @{$self->{utr5}};
}

sub get3UTR {
	my ($self) = @_;
	return @{$self->{utr3}};
}

sub cdspart {
	my ($self) = @_;
	return $self->{cds};
}

sub full {
	my ($self) = @_;
	return $self->{full};
}

sub validate {
	my ($self) = @_;
}


1;
