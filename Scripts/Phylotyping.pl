#!/usr/bin/perl
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
#
use strict;
use Bio::TreeIO;
use Getopt::Long;
use Bio::DB::Taxonomy;
$Bio::Root::Root::DEBUG = -1;

my $AMPHORA_home = $ENV{'AMPHORA2_home'};
my $tax_dir = "$AMPHORA_home/Taxonomy";
my $method = 'ml';
my (%markerlist, %taxonid) = ();
my ($CPUs, $help) = undef;
my $output = undef;

my $usage = qq~
This tool will assign each identified marker sequence a phylotype using the evolutionary placment algorithm of raxml 

Usage: $0 <options>

Options:
	-Method: use 'maximum likelihood' (ml) or 'maximum parsimony' (mp) for phylotyping. Default: ml
	-CPUs: turn on the multiple thread option and specify the number of CPUs/cores to use. If your computer has multiple CPUs/cores, you can speed up the phylotyping process by running this script on multiple cores. Important: Make sure raxmlHPC-PTHREADs is installed. If the number specified here is larger than the number of cores that are free and available, it will actually slow down the script.
	-Help: print help;  
~;


GetOptions (	'Method=s'=>\$method,
		'CPUs=i'=>\$CPUs,
		'Help'=>\$help) || die "Invalid command line options\n";

die $usage if $help;

if ($CPUs and $CPUs == 1) {
	die "CPUs has to be greater than 1";
}


my $tree_functions = new Bio::Tree::Tree(); 
my $taxdb = Bio::DB::Taxonomy->new(	-source   => 'flatfile',
					-directory=>$tax_dir,
					-nodesfile => "$tax_dir/nodes.dmp",
					-namesfile => "$tax_dir/names.dmp");

my @rank = ('species','genus','family','order','class','phylum','superkingdom');
my %rank = ();
for (@rank) {
	$rank{$_} = 1;
}

get_marker_list();
assign_phylotype();


###########################################################################################	

sub get_marker_list {
	open (IN, "$AMPHORA_home/Marker/marker.list") || die "Can't open $AMPHORA_home/Marker/marker.list";
	while (<IN>) {
		chop;
		/^(\S+)\s+(\S+)/;
		$markerlist{$1} = $2;
	}
	close IN;
}

sub assign_phylotype {
	get_contig_taxonid();
	cleanup();

	my $raxml;
	if ($CPUs) {
		$raxml = "raxmlHPC-PTHREADS -T $CPUs";
	}
	else {
		$raxml = "raxmlHPC";
	}

	
MARKER:for my $marker (keys %markerlist) {	
		next MARKER unless (-e "$marker.aln");
		if (-e "$marker.phylotype") {
			print STDERR "$marker has been assigned phylotypes; skipped...\nTo reassign phylotypes, delete the file $marker.phylotype\n\n";
			next MARKER;
		}
		my (%confidence, %support) = ();
		$output = undef;
		if ($method eq 'mp') {
			system ("$raxml -f y -t $AMPHORA_home/Marker/$marker.tre -s $marker.aln -m PROTGAMMAWAG -n $marker -p 132 1>/dev/null 2>/dev/null");
			open (IN, "RAxML_equallyParsimoniousPlacements.$marker") || die "Cannot open RAxML_equallyParsimoniousPlacements.$marker";
			while (<IN>) {
				chop;
				my ($query, $edge, $mp) = split;
				next if $query =~ /^REF-/;
				$edge =~ s/I//;
				$support{$query}{$edge} = 1;
			}
			close IN;
			
			for my $query (keys %support) {
				for my $edge (keys %{$support{$query}}) {
					$support{$query}{$edge} = 1/(scalar keys %{$support{$query}});
				}
			}	
		}
		else {
			my $count = 0;
			open (IN, "$marker.aln") ||die "cant open $marker.aln";
			while (<IN>) {
				$count++ if /REF-/; 
			}
			close IN;
			
			my $fh = 25/$count;
			$fh = 0.99 if ($fh >= 1);
			system ("$raxml -s $marker.aln -f v -G $fh -t $AMPHORA_home/Marker/$marker.tre -m PROTCATIWAG -n $marker 1>/dev/null 2>/dev/null");
			
			unless (-e "RAxML_classificationLikelihoodWeights.$marker") {
				print STDERR "Error occured when assigning phylotype for $marker\n";
				next MARKER;
			}
		
			open (IN, "RAxML_classificationLikelihoodWeights.$marker") || die "Can't open RAxML_classificationLikelihoodWeights.$marker";
			while (<IN>) {
				chop;
				my ($query, $edge, $likelihood_weight_ratio, $accumulation) = split;
				next if $query =~ /^REF-/;
				$edge =~ s/I//;
				$support{$query}{$edge} = $likelihood_weight_ratio;
			}
			close IN;
		}

		system("sed -e 's/\\[I/\\[/g'  RAxML_originalLabelledTree.$marker >  RAxML_originalLabelledTree.$marker.substitute");
		my $treein = new Bio::TreeIO('-file' =>	" RAxML_originalLabelledTree.$marker.substitute");		
		my $tree = $treein->next_tree();
		my $root = $tree->get_root_node();
		my %ref_taxon = ();
		traverse($root, \%ref_taxon);

		for my $query (keys %support) {
			$output .= "$query\t$marker";
			for my $edge (keys %{$support{$query}}) {
				my @lineage = $tree_functions->get_lineage_nodes($ref_taxon{$edge});
				for my $taxon (@lineage, $ref_taxon{$edge}) {
					$confidence{$query}{$taxon->id} += $support{$query}{$edge};	
				}
			}

			if ($markerlist{$marker} =~ /bacteria/i ) {
				assign($taxdb->get_taxon(-taxonid=>'2'), \%confidence, $query);
			}
			elsif ($markerlist{$marker} =~ /archaea/i) {
				assign($taxdb->get_taxon(-taxonid=>'2157'), \%confidence, $query);
			}

			$output .= "\n";

		}
	
		open (OUT, ">$marker.phylotype") || die "cannot open $marker.phylotype to write";
		print OUT $output;
		close OUT;
		cleanup();
	}
	print "Query\tMarker\tSuperkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
	system ("cat *.phylotype");	
}


sub get_contig_taxonid {
	open (IN, "$AMPHORA_home/Taxonomy/contig.taxonid") || die "Can't open $AMPHORA_home/Taxonomy/contig.taxonid\n";
	while (<IN>) {
		my ($contigid, $taxonid) = /^(\S+)\s+(\S+)/;
		$taxonid{$contigid} = $taxonid;
	}
	close IN;
}

sub traverse {
	my ($node, $ref_taxon) = @_;
	my (@desnodes, @taxons) = ();

	for ($node->each_Descendent) {
		push @desnodes, $_;
	}
	for (my $i=0; $i<=$#desnodes; $i++){
		my $edge = $desnodes[$i]->bootstrap();
		my $taxon;
		if (!($desnodes[$i]->is_Leaf)) {
			$taxon = traverse($desnodes[$i],$ref_taxon);
		}
		else {	
			$desnodes[$i]->id() =~ /-([^-]+)$/;
			$taxon = $taxdb->get_taxon(-taxonid=>$taxonid{$1});
		}

		$$ref_taxon{$edge} = $taxon; 
		push @taxons, $taxon;
		if ($i == $#desnodes) {
			return  $tree_functions->get_lca(\@taxons);
		}
	}
}

sub assign {
	my ($taxon,$confidence,$query) = @_;

	if ($rank{$taxon->rank}) {
		$output .= "\t".$taxon->scientific_name."(".(sprintf "%.2f", $confidence->{$query}{$taxon->id}).")";
	}
	return if ($taxon->rank =~ /species/i);
	my $best = undef;
	for my $desc ($taxdb->each_Descendent($taxon)) {
		if ($confidence->{$query}->{$desc->id}) {
			unless ($best and $confidence->{$query}->{$desc->id} < $best) {
				$best = $confidence->{$query}->{$desc->id};
				$taxon = $desc;
			}
		}
	}
	return unless $best;
	assign($taxon,$confidence,$query);
}
	 
sub cleanup {
	system ("rm RAxML* 2>/dev/null; rm *aln.reduced 2>/dev/null");
}
