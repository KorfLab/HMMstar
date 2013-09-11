#!/usr/bin/perl
use strict;
use warnings;
use HMMstar;
use Getopt::Std;
our ($opt_h, $opt_v);
getopts('hv');

die "
usage: $0 [options] <gff> <fasta>
options:
  -h  help
  -v  verbose
" unless @ARGV == 2;

my ($GFF, $DNA) = @ARGV;

my $genome = new HMMstar::Genome($DNA, $GFF);
foreach my $contig ($genome->contigs) {
	foreach my $cds ($contig->coding_sequences) {
		print $cds->name, "\t",
			substr($cds->{protein}, 0, 10), "\t",
			$cds->start_site->sequence(6, 3), "\t",
			$cds->stop_site->sequence(6, 3), "\n";
		print "\texons:";
		foreach my $exon ($cds->exons) {
			print " ", $exon->length;
		}
		print "\n";
		print "\tintrons:";
		foreach my $intron ($cds->introns) {
			print " ", $intron->length;
		}
		print "\n";
	}
}


__END__
