#! /usr/bin/env perl

use strict;
use warnings;

#declaración de variables
my $freq;
my $aafreq;
my %aaf;
my $name1;
my $codons;
my @codon_a;
my %cod;

# abre archivos de codons, <FH> en contexto escalar lee una linea. si se pone otra vez -> la siguiente
open (CODONS, "<codons.txt");
$codons = <CODONS>;
close (CODONS);
@codon_a = split(",", $codons);

if (@ARGV){
	$freq = $ARGV[0];
    $aafreq = $ARGV[1];
	open(FREQ, "<$freq");
    open(AAFREQ, "<$aafreq");
	$freq =~ /(^.+)\./;
	$freq=$1;
	$name1 = $freq.".modal_freq.fas";
	open (RESULTADOS,">$name1");
} 
else{
	print "Introduzca el nombre del archivo de frecuencias modales\n";
	$freq = <STDIN>;
    print "Introduzca el nombre del archivo de frecuencia promedio de aminoácidos\n";
    $aafreq = <STDIN>;
	open(FREQ, "<$freq");
    open(AAFREQ, "<$aafreq");
	$freq =~ /(^.+)\./;
	$freq=$1;
	$name1 = $freq.".modal_freq.fas";	
	open (RESULTADOS,">$name1");
}

chomp(my @aafile = <AAFREQ>);
(my $aarnam, my $aaname) = split("\t",$aafile[0]);
my @aanames = split(" ",$aaname);
(my $idname, my $avaafreqs) = split("\t",$aafile[1]);
my @avaafreq = split(" ",$avaafreqs);

my %aanum_h;

#hash with aa(3leters) => aafreq
my $j=0;
foreach my $aa (@aanames){
    $aanum_h{$aa} = $avaafreq[$j]*10000;
    $j++;
}

#declarar variables de todo el namespace-package
my %aaorder=(
0 => "Ala",
1 => "Ala",
2 => "Ala",
3 => "Ala",
4 => "Cys",
5 => "Cys",
6 => "Asp",
7 => "Asp",
8 => "Glu",
9 => "Glu",
10 => "Glu",
11 => "Glu",
12 => "Phe",
13 => "Phe",
14 => "Gly",
15 => "Gly",
16 => "His",
17 => "His",
18 => "Ile",
19 => "Ile",
20 => "Ile",
21 => "Lys", 
22 => "Lys",
23 => "Leu",
24 => "Leu",
25 => "Leu",
26 => "Leu",
27 => "Leu",
28 => "Leu",
29 => "Asn",
30 => "Asn",
31 => "Pro",
32 => "Pro",
33 => "Pro",
34 => "Pro",
35 => "Gln",
36 => "Gln",
37 => "Arg",
38 => "Arg",
39 => "Arg",
40 => "Arg",
41 => "Arg",
42 => "Arg",
43 => "Ser",
44 => "Ser",
45 => "Ser",
46 => "Ser",
47 => "Ser",
48 => "Ser",
49 => "Thr",
50 => "Thr",
51 => "Thr",
52 => "Thr",
53 => "Val",
54 => "Val",
55 => "Val",
56 => "Val",
57 => "Tyr",
58 => "Tyr",
59 => "Met",
60 => "Trp");

while (<FREQ>) {
	my $line = $_;
	chomp $line;
	$line =~ s/^>//;
	$line =~ s/\|/,/g;
	(my $lineid, my $modal) = split("\t",$line);
	my @modal_a = split(",",$modal);
	
    my $i=0;
    my $max_index = $#modal_a;
    if ($max_index < 60){
        push @modal_a, 1;
        push @modal_a, 1;
    }
	foreach my $x (@modal_a) {
        $x = $x * $aanum_h{$aaorder{$i}}; 
        $i++;
	}	
    $max_index = $#modal_a;
	for (my $i=0; $i <= $max_index; $i++){
		$cod{$codon_a[$i]} = $modal_a[$i];
	}
    

	#algoritmo
	print RESULTADOS ">".$freq."\n"."ATG";
	foreach my $key (keys %cod){
    	print RESULTADOS "$key"x(int($cod{$key}+0.5));
	}
	print RESULTADOS "TGA\n";
}
close FREQ;
close AAFREQ;
close RESULTADOS;

