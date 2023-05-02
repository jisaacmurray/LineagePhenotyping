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
    my $cmd = "cp " . $DBref->{$series}{"annots"} . "/dats/*.csv $target";
    system($cmd);
}

close(FILE)
