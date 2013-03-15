# AMPHORA (version 2.0) An Automated Phylogenomic Inference Pipeline for Bacterial and Archaeal Sequences. 
# Copyright 2011 by Martin Wu
 
# This file is part of AMPHORA2.

# AMPHORA2 is free software: you may redistribute it and/or modify its under the terms of the 
# GNU General Public License as published by the Free Software Foundation; either version 2 of
# the License, or any later version.

# AMPHORA2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details (http://www.gnu.org/licenses/).
 

# For any other inquiries send an Email to Martin Wu
#       mw4yv@virginia.edu
 
# When publishing work that is based on the results from AMPHORA2 please cite:
# Wu M, Scott AJ: Phylogenomic analysis of bacterial and archaeal sequences with AMPHORA2. Bioinformatics 2012, 28:1033-1034.

use strict;
use Bio::SeqIO;
use Getopt::Long;

my $AMPHORA_home = $ENV{'AMPHORA2_home'};
my (%markerlist, %seq) = ();
my $help = undef;
my $evalue = 1e-7;
my $ref_dir = "$AMPHORA_home/Marker/";
my ($is_dna, $is_bacteria, $is_archaea) = undef;
$Bio::Root::Root::DEBUG = -1;

my $usage = qq~
This tool will search for bacterial and archaeal phylogenetic markers for a given fasta sequence file.
The output of this tool is a collection of marker protein sequences in fasta format. For example, rpoB.pep, dnaG.pep.

When DNA sequences are used, this program first identifies ORFs longer than 100 bp in all six reading frames, then scans the translated peptide sequences for the phylogenetic markers.

The default is assuming that the query sequences contain both bacterial and archaeal sequences. If you know your query sequences only contain bacterial sequences, then use the -Bacteria option. If you know your query sequences only contain archaeal sequences, then use the -Archaea option. This makes the program run faster and the results will be more accurate.

Usage: $0 <options> sequence-file

Options:
	-DNA: input sequences are DNA. Default: no.
	-Evalue: HMMER evalue cutoff. Default: 1e-7 
	-Bacteria: input sequences are Bacterial sequences
	-Archaea: input sequences are Archaeal sequences
	-ReferenceDirectory: the file directory that contain the reference alignments, hmms and masks. Default: $AMPHORA_home/Marker
	-Help: print help;
~;



GetOptions (	'DNA'=>\$is_dna,
		'Evalue=f'=>\$evalue,
		'Bacteria'=>\$is_bacteria,
		'Archaea'=>\$is_archaea,
		'ReferenceDirectory=s'=>\$ref_dir,
		'Help'=>\$help) || die "Invalid command line options\n";

die $usage if $help;
die $usage unless $ARGV[0];

get_marker_list();

if ($is_bacteria) {
	$ref_dir .= "Bacteria.";
}
elsif ($is_archaea) {
	$ref_dir .= "Archaea.";
}

my $input_seq = $ARGV[0];

if ($is_dna) {
	system ("getorf -sequence $ARGV[0] -outseq $ARGV[0].orf -table 1 -minsize 100");
	$input_seq = "$ARGV[0].orf";
}

# HMM search
system ("hmmsearch -Z 5000 -E $evalue --domE $evalue --domtblout $$.hmmsearch -o /dev/null $ref_dir"."markers.hmm $input_seq");		# fix the number of sequences in the database for E-value calculation
get_hmm_hits();

# clean up
system ("rm $$.*");

####################################################################################################################
sub get_marker_list {
	open (IN, "$AMPHORA_home/Marker/marker.list") || die "Can't open $AMPHORA_home/Marker/marker.list";
	while (<IN>) {
		chop;
		/^(\S+)\s+(\S+)/;
		$markerlist{$1} = $2;
	}
	close IN;
}

sub get_hmm_hits {
	my (%score, %candidates, %query_length, %hmm_length, %query_match, %hmm_match, %hits) =();
	open (IN, "$$.hmmsearch") || die "Can't open $$.hmmsearch";
	while (<IN>) {
		chop;
		next if /^#/;
		my ($query, $query_accession, $qlength, $hmm, $hmm_accession, $hmm_length, $evalue, $score, $bias, $domain, $domain_number, $dom_evalue, $ievalue, $dom_score, $dom_bias, $hmm_start, $hmm_stop, $query_start, $query_stop, $rest) = split /\s+/;
		if ((! exists $score{$query}) or ($score{$query} < $score)) {
			$score{$query} = $score;
			$candidates{$query} = $hmm;
		}
		$query_length{$query} = $qlength;
		$hmm_length{$hmm} = $hmm_length;
		$query_match{$hmm}{$query} += ($query_stop - $query_start);
		$hmm_match{$hmm}{$query} += ($hmm_stop - $hmm_start);
	}
	close IN;
	
	while ( my ($query, $hmm) = each (%candidates) ) {
		next unless $markerlist{$hmm};
		next if ( ($query_match{$hmm}{$query}/$query_length{$query} < 0.7) and ($hmm_match{$hmm}{$query}/$hmm_length{$hmm} < 0.7) ); 	#ignore the hit if the match is partial
		$hits{$hmm}{$query} = 1;
	}

	my $seqin = new Bio::SeqIO('-file'=>$input_seq);
	while (my $seq = $seqin->next_seq) {
		$seq{$seq->id} = $seq;
	}
	
	for my $marker (keys %hits) {
		my $seqout = new Bio::SeqIO('-file'=>">>$marker.pep",'-format'=>'fasta');
		for my $seqid (keys %{$hits{$marker}}) {
			$seqout->write_seq($seq{$seqid});
		}
	}
}
