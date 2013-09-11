package HMMstar::Genome;
use strict;
use warnings 'FATAL' => 'all';
use HMMstar::CDS;
use HMMstar::CodingTranscript;
use HMMstar::Contig;
use HMMstar::Decoder;
use HMMstar::Feature;
use HMMstar::HMM;
use HMMstar::SequenceModel;
use HMMstar::State;
use HMMstar::Transcript;

my $DEBUG = 0;

sub new {
	my ($class, $fasta, $gff) = @_;
	
	my $sequences = read_fasta($fasta);	
	my $features  = read_gff($gff, $sequences);
	
	my $self = bless {
		sequences => $sequences,
		features  => $features,
	};
	$self->validate;
	return $self;
}

sub validate {
	my ($self) = @_;
		
	# sanity check: make sure sequence identifiers in gff are in fasta
	my @unique_keys;
	foreach my $key (keys %{$self->{features}}) {
		if (not exists $self->{sequences}{$key}) {
			push @unique_keys, $key;
		}
	}
	if (@unique_keys) {
		die "Some sequence identifiers in GFF not in FASTA: @unique_keys\n";
	}
}

sub read_gff {
	my ($file, $seqs) = @_;
		
	my %gff;
	open(IN, $file) or die;
	while (<IN>) {
		next unless /\S/;
		next if /^#/;
		my @f = split;
		push @{$gff{$f[0]}}, new HMMstar::Feature($seqs->{$f[0]},
			$f[2], $f[3], $f[4], $f[6], $f[5], $f[7]);
	}
	close IN;
	
	if ($DEBUG) {
		print STDERR "read ", scalar keys %gff, " annotations\n";
	}
	
	return \%gff;
}

sub read_fasta {
	my ($file) = @_;
		
	my %seq;
	my $id;

	
	open(IN, $file) or die;
	while (<IN>) {
		next unless /\S/;
		if (/^>(\S+)/) {$id = $1}
		else           {chomp; $seq{$id} .= $_}
	}
	
	foreach my $id (keys %seq) {
		$seq{$id} = new HMMstar::DNA($id, $seq{$id});
	}
	
	if ($DEBUG) {
		print STDERR "read ", scalar keys %seq, " sequences\n";
	}
	
	return \%seq;
}

sub contigs {
	my ($self) = @_;
	
	my @contig;
	foreach my $id (sort keys %{$self->{sequences}}) {
		push @contig, new HMMstar::Contig($self->{sequences}{$id},
			$self->{features}{$id});
	}
	return @contig;
}

1;
