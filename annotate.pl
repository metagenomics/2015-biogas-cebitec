#!/usr/bin/env perl

# Some random Perl script to combine BLASTP and BEDTools output with KEGG mappings
# Requires the KEGG mapping files 'genes_ko.list', 'ko_module.list', and 'ko_pathway.list'
# These, together with BLASTP and BEDTools output files, must be in the working directory
# For a sample output, please see: Contigs_gt1kb.prodigal.faa.blastp.kegg.annotated.tsv.gz
#
# Author: abremges@cebitec.uni-bielefeld.de

use strict;
use warnings;

my $file;
my %hash = ();
my %genes_ko = ();
my %ko_module = ();
my %ko_pathway = ();

open($file, '<:encoding(UTF-8)', 'genes_ko.list') or die $!;
while (my $line = <$file>) {
	chomp $line;
	my @tmp = split('\t', $line);
	push @{ $genes_ko{$tmp[0]} }, $tmp[1];
}
close $file;

open($file, '<:encoding(UTF-8)', 'ko_module.list') or die $!;
while (my $line = <$file>) {
	chomp $line;
	my @tmp = split('\t', $line);
	push @{ $ko_module{$tmp[0]} }, $tmp[1];
}
close $file;

open($file, '<:encoding(UTF-8)', 'ko_pathway.list') or die $!;
while (my $line = <$file>) {
	chomp $line;
	my @tmp = split('\t', $line);
	push @{ $ko_pathway{$tmp[0]} }, $tmp[1];
}
close $file;

open($file, '<:encoding(UTF-8)', '4_genes/Contigs.prodigal.gff.bedtools.tsv') or die $!;
while (my $line = <$file>) {
	chomp $line;
	$line =~ m/^(contig-\d+).*ID=(\d+)_(\d+)/;
	my $gene_id = $1 . "_" . $3;
	my @tmp = split('\t', $line);
	$hash{$gene_id}{'len'} = $tmp[4]-$tmp[3];
	$hash{$gene_id}{'dna'} = $tmp[10]+$tmp[11]+$tmp[12]+$tmp[13]+$tmp[14]+$tmp[15];
	$hash{$gene_id}{'rna'} = $tmp[9];
}
close $file;

open($file, '<:encoding(UTF-8)', '4_genes/Contigs.prodigal.faa.blastp.kegg.tsv') or die $!;
while (my $line = <$file>) {
	chomp $line;
	my @tmp = split('\t', $line);
	$hash{$tmp[0]}{'gene'} = $tmp[1];
	$hash{$tmp[0]}{'eval'} = 0+$tmp[10];
}
close $file;

print "gene_id\tgene_length\tdna_reads\trna_reads\tkegg_gene\tblastp_evalue\tkegg_ko\tkegg_module\tkegg_pathway\n";
foreach my $key ( sort keys %hash ) {
	print "$key\t";

	if (exists($hash{$key}{'len'})) {
		print "$hash{$key}{'len'}\t";
	} else {
		print "*\t";
	}

	if (exists($hash{$key}{'dna'})) {
		print "$hash{$key}{'dna'}\t";
	} else {
		print "*\t";
	}

	if (exists($hash{$key}{'rna'})) {
		print "$hash{$key}{'rna'}\t";
	} else {
		print "*\t";
	}

	if (exists($hash{$key}{'gene'})) {
		print "$hash{$key}{'gene'}\t";
	} else {
		print "*\t";
	}

	if (exists($hash{$key}{'eval'})) {
		print "$hash{$key}{'eval'}\t";
	} else {
		print "*\t";
	}

	if (exists($hash{$key}{'gene'}) && exists($genes_ko{$hash{$key}{'gene'}})) {
		print join(",", sort @{$genes_ko{$hash{$key}{'gene'}}});
		print "\t";
	} else {
		print "*\t";
	}

	if (exists($hash{$key}{'gene'}) && exists($genes_ko{$hash{$key}{'gene'}})) {
		my %tmp = ();
		foreach my $ko (@{$genes_ko{$hash{$key}{'gene'}}}) {
			foreach my $module (@{$ko_module{$ko}}) {
				$tmp{$module} = 1;
			}
		}
		if (keys %tmp) {
			print join(",", sort keys %tmp);
			print "\t";
		} else {
			print "*\t";
		}
	} else {
		print "*\t";
	}

	if (exists($hash{$key}{'gene'}) && exists($genes_ko{$hash{$key}{'gene'}})) {
		my %tmp = ();
		foreach my $ko (@{$genes_ko{$hash{$key}{'gene'}}}) {
			foreach my $pathway (@{$ko_pathway{$ko}}) {
				$tmp{$pathway} = 1;
			}
		}
		if (keys %tmp) {
			print join(",", sort keys %tmp);
			print "\n";
		} else {
			print "*\n";
		}
	} else {
		print "*\n";
	}
}
