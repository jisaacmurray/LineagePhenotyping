#!/usr/bin/perl -w
use strict;
use MakeDB;

my $DBref = MakeDB();
my $list = $ARGV[0];
my $target = $ARGV[1];
#argument: list file, target directory
open(FILE, $list);
my $series;
foreach $series (<FILE>){
    chomp $series;
    my $CD = $DBref->{$series}{"annots"} . "/dats/CD" . $series . ".csv";
    print $CD . "\n";

    if(-e $CD){
	foreach(`cat $CD`){
	    print "$series:$_";
	}
    }else{
	print STDERR "NO FILE $CD\n";
    }
    
}

close(FILE)
