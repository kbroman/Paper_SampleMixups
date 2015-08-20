#!/usr/bin/perl
$ifile = "OrigData/Final Fit1_Filtered_Assay_Allele_Signals_and_Genotypes_18Sep.txt";
$ofile = "R/snp_intensities.csv";
$ofile2 = "R/snp_averages.csv";

open(IN, $ifile) or die("Cannot read from $ifile");
open(OUT, ">$ofile") or die("Cannot write to $ofile");
print OUT "MouseNum,A,C,G,T,SignalSum\n";
open(OUT2, ">$ofile2") or die("Cannot write to $ofile2");
print OUT2 "MouseNum,AveSignalSum\n";

$line = <IN>;
while($line = <IN>) {
    chomp($line);
    @v = split(/\t/, $line);
    $id = $v[4];
    if($id =~ /_Mouse(\d\d\d\d)_/) {
        $id = "Mouse$1";
    }
    @sig = @v[13..16];
    $sigsum = $v[18];
    print OUT (join(",", ($id,@sig,$sigsum)), "\n");
}
