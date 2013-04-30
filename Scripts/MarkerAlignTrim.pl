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
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Getopt::Long;
$Bio::Root::Root::DEBUG = -1;


my $AMPHORA_home = $ENV{'AMPHORA2_home'};

my $usage = qq~
Usage:  $0 <options>

Assume the Phylogenetic Marker Database is located at $AMPHORA_home

Options:
	-Trim:		trim the alignment using masks embedded with the marker database. Default: no
	-Cutoff:	the Zorro masking confidence cutoff value (0 - 1.0; default: 0.4);
	-ReferenceDirectory: the file directory that contain the reference alignments, hmms and masks. Default: $AMPHORA_home/Marker
	-Directory:	the file directory where sequences to be aligned are located. Default: current directory
	-OutputFormat:  output alignment format. Default: phylip. Other supported formats include: fasta, stockholm, selex, clustal
	-WithReference:  keep the reference sequences in the alignment. Default: no
	-Help:		print help 
~;

my ($trim, $help, $with_ref) = undef;
my $dir = ".";
my $cutoff = 0.4;
my $ref_dir = "$AMPHORA_home/Marker";
my $format = 'phylip';
my $alignment;

my (%markerlist, %imprint, %mask, @output_seq) = ();

GetOptions (	'Trim'=>\$trim,
		'Cutoff=f'=>\$cutoff,
		'ReferenceDirectory=s'=>\$ref_dir,
		'Directory=s'=>\$dir,
		'OutputFormat=s'=>\$format,
		'WithReference'=>\$with_ref,
		'Help'=>\$help) || die "Invalid command line options\n";

die $usage if $help;
die $usage unless ("-e $ref_dir/marker.list");

get_marker_list();

for my $marker (keys %markerlist) {
	next unless (-e "$dir/$marker.pep");
	print STDERR "Aligning $marker ...\n";
	$alignment = new Bio::SimpleAlign();
	(%imprint, %mask, @output_seq) = ();
	read_mask($marker);
	align($marker);
	trim();
	output($marker);
	clean();
}



####################################################################################################
sub get_marker_list {
	open (IN, "$ref_dir/marker.list") || die "Can't open $ref_dir/marker.list";
	while (<IN>) {
		chop;
		/^(\S+)/;
		$markerlist{$1} = 1;
	}
	close IN;
}
			
sub read_mask {
	my $marker = shift;
	my @weight = ();
	open (IN, "$ref_dir/$marker.mask") || die "Can't open $ref_dir/$marker.mask";
	while (<IN>) {
		chop;
		push @weight, $_;
	}
	close IN;

	open (IN, "$ref_dir/$marker.stock") || die "Can not open $ref_dir/$marker.stock";
	while (<IN>) {
		chop;
		next if /^#/;
		my $j = 0;
		my ($id, $seq) = /^(\S+)\s+(\S+)$/;
		for (my $i=0; $i< length($seq); $i++) {
		next if (substr($seq, $i, 1) eq '-');
            		if ($i <= $#weight) {
				$imprint{$id}{$j} = $weight[$i]/10;
		    	}
		    	else {
				$imprint{$id}{$j} = '0';
			}
			$j++;
		}
	}	

	my $seqin = new Bio::SeqIO('-file'=>"$dir/$marker.pep");
	while (my $seq = $seqin->next_seq) {
		die "Query sequence ".$seq->id()." is present in the reference sequences\n" if exists $imprint{$seq->id()}; 
	}
}

	
sub align {
	my $marker = shift;

	# remove duplicated sequences in the alignment
	remove_duplicate($marker);

	system("hmmalign -o $$.slx --mapali $ref_dir/$marker.stock $ref_dir/$marker.hmm $dir/$marker.pep>/dev/null"); 	
	my $in = new Bio::AlignIO('-file'=>"$$.slx");
	my $alignment = $in->next_aln();
	my @seq = $alignment->each_seq();
	
	for my $seq (@seq) {
		my ($ID) = ($seq->id() =~ /^(\S+)/);
		$seq->id($ID);
		my $sequence = $seq->seq();
		$sequence =~ s/\./-/g;
		$seq->seq($sequence);
		if (exists $imprint{$ID}) {
			my $j=0;
			LOOP: for (my $i=0; $i<$seq->length(); $i++) {
				next LOOP if (substr($sequence, $i, 1) eq '-');
				$mask{$i} = $imprint{$ID}{$j} if (exists $imprint{$ID}{$j});
				$j++;
			}
			if ($with_ref) {
				push @output_seq, $seq;
			}
		}
		else {
			push @output_seq, $seq;
		}
	}
}

sub trim {
	for my $seq (@output_seq) {
		if ($trim) {
			my $i = 0;
			my $sequence = '';
			for (split //,$seq->seq()) {
				$sequence .= $_ if ( (exists $mask{$i}) and ($mask{$i} >= $cutoff) );
				$i++;
			}
			$seq->seq($sequence);
		}
		$alignment->add_seq($seq);
	}
}	

sub remove_duplicate {
	my $marker = shift;
	my %unique =();
	my $seqin = new Bio::SeqIO('-file'=>"$dir/$marker.pep");
	while (my $seq = $seqin->next_seq) {
		$unique{$seq->id()} = $seq;
	}
	my $seqout = new Bio::SeqIO('-file'=>">$dir/$marker.pep", '-format'=>'fasta');
	for (keys %unique) {
		$seqout->write_seq($unique{$_});
	}
}

sub output {
	my $marker = shift;
	my $alignout;
	if ($format =~ /phylip/i) {
		$alignout = new Bio::AlignIO('-file'=>">$dir/$marker.aln",'-format'=>$format, '-idlength'=> 50);
	}
	else {
		$alignout = new Bio::AlignIO('-file'=>">$dir/$marker.aln",'-format'=>$format);
	}
	$alignout->write_aln($alignment);

	open (OUT, ">$dir/$marker.mask") || die "Can't write $dir/$marker.mask\n";
	my @tmp = sort {$b <=> $a} keys %mask;
	for (1 .. $tmp[0]) {
		if (exists $mask{$_}) {
			print OUT int($mask{$_}*10)," ";
		}
		else {
			print OUT "0 ";
		}
	}
	print "\n";
	close OUT;
}

sub clean {
	system ("rm $$.*");
}
