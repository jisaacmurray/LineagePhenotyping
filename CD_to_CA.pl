#!/usr/bin/perl -w
use strict;

##this is specifically for Ma et al CD files from /murrlab3/jmurr/Ma_Du/CD


print join(",", qw(celltime cell time none global local blot cross z x y size gweight)) . "\n";
my $ct=0;
my $exp=0;
my $prevCell = "NULL";
my @prevData;
foreach(<>){
    chomp;
#    print;
    my @data = split(/,/);
    my $cell = $data[1];
    my $e1 = $data[3];
#    print "$cell $e1\n";
    if($cell eq $prevCell){
	$ct = $ct+1;
	$exp = $exp+$e1
    }else{
	if($prevCell ne "NULL"){
	    my $total = int($exp/$ct*10);
	    print join(",", @prevData[0..2], $total, $total, $total, $total, $total, @prevData[8..12]) . "\n";
	    $ct=0;
	    $exp=0;
	}
    }
    $prevCell = $cell;
    @prevData = @data;
}

my $total = int($exp/$ct*10);
print join(",", @prevData[0..2], $total, $total, $total, $total, $total, @prevData[8..12]) . "\n";
