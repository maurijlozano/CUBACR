#!/usr/bin/perl -w
#Usage: perl codonAlnStats.pl

#Citation: Lozano M.J
#V1. 2019

use strict;

sub usage {
	my $error = shift;
	print STDERR " codonAlnStats.pl (2019)\n\n";
	print STDERR " Error: $error\n\n";
	print STDERR "             Lozano Mauricio J.\n";
	print STDERR "             Transforms Multiople alignment to table .\n";
	print STDERR "             and extracts conserved an variable positions.\n\n\n";
	print STDERR " Basic usage:    codonAlnStats.pl -i nt_file\n";
	print STDERR "  -i: the file containing the nucleotide alignment, fasta format. \n";
	print STDERR "  -o: output file (Optional).\n";
	print STDERR "  -a: file containing the amino acid sequence alignment (Optional)\n";
	print STDERR "\n";
	exit;
}

my %params = @ARGV;
my %valid_params;
$valid_params{"-i"} = 1;	$valid_params{"-o"} = 1;	$valid_params{"-a"} = 1;	

foreach my $param ( keys %params ) {
	if(!exists($valid_params{$param})) {
		&usage("Option -$param not recognized");
	}
	if(!defined($params{$param})) {
		&usage("Option -$param requires an accompanying value");
	}
}


if(!exists($params{"-i"})) {	&usage("Option -i not defined") ;	}
my $nt_file = $params{"-i"};
if(!-e $nt_file) {	&usage("File $nt_file does not exist");		}

if(!exists($params{"-a"})) {	&usage("Option -a not defined") ;	}
my $aa_file = $params{"-a"};
if(!-e $aa_file) {	&usage("File $aa_file does not exist");		}


my $out_file = "";
if(exists($params{"-o"})) {	$out_file = $params{"-o"};	}
else{	$out_file = "ali2table";		}


# Variables.
my @nt_Seq_id;
my @nt_Seq_entry;
my @nt_seq;
my @aa_Seq_id;
my @aa_Seq_entry;
my @aa_seq;
my $aa_out_file;
my $nt_out_file;
my $id;

# #########################################################################
# 
# Alignment to table
# 
# 
# #########################################################################


# Read in nucleotide alignment and transform to table
open (FILE1, "<$nt_file")  or die "Could not open file '$nt_file'. $!";
while (<FILE1>) {
	chomp;
	if ($_ =~ /^>/){
		$id = $_;
		$id =~ s/^>//;
		push @nt_Seq_id, $id;
	}
	else{
		push @nt_Seq_entry, $_;
	}
}	
close FILE1;	

open (FILE2, "<$aa_file")  or die "Could not open file '$aa_file'. $!";
while (<FILE2>) {
	chomp;
	if ($_ =~ /^>/){
		$id = $_;
		$id =~ s/^>//;
		push @aa_Seq_id, $id;
	}
	else{
		push @aa_Seq_entry, $_;
	}
}	

close FILE2;	

my $len   = @nt_Seq_entry;
# Split aa or codons and save table
$nt_out_file = $out_file.".nt.table";
open (FILE3, ">$nt_out_file");
for (my $i=0; $i < $len; $i++){
	@nt_seq = $nt_Seq_entry[$i] =~ /(...)/g;
	print FILE3 $nt_Seq_id[$i]." ".join(" ", @nt_seq)."\n";
}
close FILE3;


$len = @aa_Seq_entry;
$nt_out_file = $out_file.".aa.table";
open (FILE4, ">$nt_out_file");
for (my $i=0; $i < $len; $i++){
	@aa_seq = $aa_Seq_entry[$i] =~ /(.)/g;
	print FILE4 $aa_Seq_id[$i]." ".join(" ", @aa_seq)."\n";
}
close FILE4;



