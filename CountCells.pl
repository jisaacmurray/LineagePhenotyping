#!/usr/bin/perl
use strict;
use MakeDB;


my $DBref = MakeDB();
my @series;
foreach(<>){
    chomp;
    push(@series, $_);
}

my $series;
print "series\tcount\tfraction\n";
foreach $series(@series){
    my $dir = $DBref->{$series}{"annots"} . "/dats/";
    my $CA = $dir . "CA" . $series . ".csv";
    if(-e $CA){
	my $count = 0;
	foreach(`cat $CA`){
	    $count++;
	}
	print "$series\t$count\t" . $count/1341 . "\n";
    }else{
	print STDERR "Can't find $CA\n";
    }
}
