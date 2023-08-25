#!/usr/bin/perl

use warnings;
use strict;

if (@ARGV != 2) {
	print "usage: $0 <input> <output prefix>\n";
	print "where:\n\t<input> name for the input file\n";
	print "\t<output> prefix for the per locus output fasta files\n";
	exit;
}

chomp (my $inputFile = $ARGV[0]);
chomp (my $outputPrefix = $ARGV[1]);

open (IN, $inputFile), or die "There was a problem opening the input file: $inputFile\n$!\n";

my $hdr = <IN>;
chomp $hdr;

my @colList = split(/\t/, $hdr, -1);
my %cols = ();
for (my $i = 0; $i < scalar @colList; $i++) {
	my $cols = $colList[$i];

	$cols{$cols} = $i;
}

my %seqs = ();

while (my $l = <IN>) {
	chomp $l;

	my @data = split(/\t/, $l, -1);

	my $locus = $data[($cols{locus})];
	my $id = $data[($cols{sequence_id})];
	my $seq = $data[($cols{sequence})];


	$seqs{$locus}{$id} = $seq;
}
close IN;

#generate the output file for each locus
foreach my $locus (keys %seqs) {
	my $outputFile = $outputPrefix . "_" . $locus . ".fa";
	open (OUT, ">$outputFile"), or die "There was a problem opening the output file: $outputFile\n$!\n";

	foreach my $id (keys %{$seqs{$locus}}) {
		my $seq = $seqs{$locus}{$id};

		print OUT ">" . $id . "\n" . $seq . "\n";
	}

	close OUT;
}

#done
exit;