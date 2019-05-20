#! /usr/bin/env perl
# Programed by Mauricio J Lozano
use strict;
use warnings;

#variables
my %cod;
my %aa;
my @codons_a;
my @seqcodons;
my $tempfile;
my $counts;
my $head;
my $head2;




# abre archivos de barcodes y pools
open (CODONS, "<codons2.txt");
@codons_a = <CODONS>;
close (CODONS);

#indexa los codones, el número de cada codones por aa etc
foreach (@codons_a) {
	(my $aas, my $ncod, my $codon) = split(",", $_);
	chomp($codon);
	$cod{$aas}=$codon;
    $aas =~ /\d\d_(...)_.../;
	$aa{$1}=$1;
}

$head = "Genes/Codons";
foreach (sort keys %cod){
	$head = $head." ".$_;
}

$head2 = "AA\t";
foreach (sort keys %aa){
	$head2 = $head2.$aa{$_}." ";
}
$head2 =~ s/ $//;

#archivos
my $file1;

if (@ARGV){
	$file1 = $ARGV[0];
	} 
else{
    print "Please type in the file name: ";
    $file1 = <STDIN>;
}

$file1 =~ /(^.+)\./;
my $name = $1;
print $name;
my $file2 = $name.".counts";
my $file3 = $name.".aaav";

#input.dat and temp.file
$tempfile = readfile ($file1,"seq.temp");

#calculate counts
$counts = counts ($tempfile,$file2,$file3);


unlink "seq.temp";


#**********************************************
sub readfile {
	my @files=@_;
	
	my $input=$files[0];
    my $temp=$files[1];
	my @tempo;
	my $i=0;

	open (TEMP,">$temp");
	open (DATOS, "<$input");
	chomp(my @lines = <DATOS>);

	 foreach my $line (@lines) {
		if ($line =~ /^>/){
			$line =~ s/\s/_/g;
			$line =~ s/\'/_/g;
			push @tempo, "\n".$line."\n";
		}
		elsif($line =~ /^$/){
		}
		else {
			$line =~ s/\r//g;
            $line =~ s/-*$//;
			$line =~ tr/agtc/AGTC/;
			push @tempo, $line;
			}
	}
	$tempo[0] =~ s/^\n//;
	foreach my $line (@tempo){
		print TEMP $line;
	}
	close TEMP;
	close DATOS;
	return $temp;
}


sub counts {
	my @files=@_;
	
	my $input=$files[0];
    my $count=$files[1];
    my $aaav=$files[2];
	my $i=0;
	my $id;
	my $nseqs=0;
	my %aafreqsum;
	my %aaav;
	my %res_counts;

	open (AAAV,">$aaav");
	open (DATOS, "<$input");

	$res_counts{0} = $head;
    print AAAV $head2."\n";
    print AAAV $name."\t";
    chomp(my @lines = <DATOS>);


    #calculate counts
	foreach my $line (@lines){
		my %codoncount;
		my %codoncount2;
		my @seqcount;

		if ($line =~ /^>/){
			$line =~ s/>//;
			$id = $line;
	        $i++;
		}
		else {
			foreach my $acod (sort keys %cod){
				$codoncount{$acod} = 0;
				$codoncount2{$cod{$acod}}=0;
			}

			@seqcodons = $line =~ /(...)/g;
			my $codons_len = @seqcodons;
			if ($codons_len >2){
				pop @seqcodons;
		        shift @seqcodons;	
            }
			foreach my $codon (@seqcodons){
				$codoncount2{$codon}++;
			}
			
			foreach (sort keys %cod){
				my $acod=$_;
				$codoncount{$acod} = $codoncount2{$cod{$acod}};
			}

			foreach my $keys (sort keys %codoncount){
				push @seqcount,$codoncount{$keys};
			}

			my $counted = join(' ',@seqcount);
			$res_counts{$i} = $id." ".$counted;

			#calculate aaAV
			$nseqs++;

   			my @seqaaav;
   			my %aacount;
            my $aalen;
			my %codfreqs = %codoncount;

			foreach my $key (sort keys %codfreqs){
				$key =~ /.{3}(\w\w\w)/;
				my $key1 = $1;
				$aacount{$key1} += $codfreqs{$key};
                $aalen += $codfreqs{$key};
			}

	    	foreach my $key (sort keys %aacount){
                $aafreqsum{$key}+=$aacount{$key}/$aalen;
			}
		}
	}
	my $printres="";

    foreach my $key (sort keys %aafreqsum){
				$aaav{$key} = $aafreqsum{$key}/$nseqs;
                $printres = $printres.$aaav{$key}." ";
  	}
    $printres =~ s/ $//;
    print AAAV $printres;

	return 0;

	close AAAV;
	close DATOS;
}




























