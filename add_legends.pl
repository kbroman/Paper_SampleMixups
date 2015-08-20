#!/usr/bin/perl
#
# add legends to files

$ifile = "legends.txt";
open(IN, $ifile) or die("Cannot read from $ifile");
while($line = <IN>) {
    chomp($line);
    if($line =~ /^>(KWB_[A-Za-z0-9]*_LEGEND)/) {
        $tag = $1;
        $text = "";
        while($line = <IN>) {
            if($line =~ /^</) { last; }
            $text .= $line;
        }
        $legend{$tag} = $text;
    }
}

$ifile = $ARGV[0];
$ofile = $ARGV[1];
open(IN, $ifile) or die("Cannot read from $ifile");
open(OUT, ">$ofile") or die("Cannot write to $ofile");
while($line = <IN>) {
    if($line =~ /(KWB_[A-Za-z0-9]+_LEGEND)/) {
        $tag = $1;
        $line =~ s/$tag/$legend{$tag}/;
    }
    print OUT $line;
}
