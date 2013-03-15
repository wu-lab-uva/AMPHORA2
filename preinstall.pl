#! /usr/bin/perl -w
use strict;
use LWP::Simple;

my $path = `pwd`;

#Check if HMMER3 is installed
my $hmm_config = `hmmsearch -h`;

if($hmm_config =~ /HMMER 3\.0/){
	print "HMMER3 is installed.\n";
}
else{
	print "HMMER3 not found, installing HMMER3 now.\n";
	#Install HMMER3
	getstore('http://selab.janelia.org/software/hmmer3/3.0/hmmer-3.0.tar.gz', "hmmer-3.0.tar.gz");
	system("tar xf hmmer-3.0.tar.gz");
	chdir "hmmer-3.0";
	system("./configure --prefix=/usr/local/");
	system("make");
	system("make install");
}
chdir;

#Check if RAxML installed
my $raxml_config = `which raxmlHPC`;
my $raxml_PTHREADS_config = `which raxmlHPC-PTHREADS`;

if(($raxml_config =~ /.+raxmlHPC/) and ($raxml_PTHREADS_config =~ /.+raxmlHPC-PTHREADS/) ){
	print "RAxML is installed.\n";
}
else{
	print "RAxML not found, installing RAxML now.\n";
	#Install RAxML
	getstore('https://github.com/stamatak/standard-RAxML/tarball/master', "stamatak-standard-RAxML.tar.gz");
	system("tar xf stamatak-standard-RAxML.tar.gz");
	system ("cd stamatak-standard-RAxML-*; make -f Makefile.gcc; make -f Makefile.PTHREADS.gcc; cp raxmlHPC* /usr/local/bin/.");
}
chdir;

no warnings 'uninitialized';
#Check if BioPerl is installed
my $perl_config = `perl -MBio::Root::Version -e 'print \$Bio::Root::Version::VERSION'`;
if($perl_config =~ /\d.+/){
	print "BioPerl is installed.\n";
}
else{
	print("BioPerl not found, installing BioPerl now.\n");
	#Install BioPerl
	getstore('http://bioperl.org/DIST/BioPerl-1.6.1.tar.gz', "BioPerl-1.6.1.tar.gz");
	system("tar xvfz BioPerl-1.6.1.tar.gz");
	chdir "BioPerl-1.6.1";
	system("perl Build.PL");
	system("./Build install");
}
chdir;

#Check if EMBOSS is installed
my $emboss_config = `which getorf`;
if($emboss_config =~ /.+getorf/){
	print "getorf is installed.\n";
}
else{
	print "getorf not found, installing EMBOSS now.\n";
	#Install EMBOSS
	getstore('ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.5.7.tar.gz', "EMBOSS-6.5.7.tar.gz");
	system("tar xf EMBOSS-6.5.7.tar.gz");
	chdir "EMBOSS-6.5.7";
	system("./configure");
	system("make");
	system("make install");
}
chdir;
