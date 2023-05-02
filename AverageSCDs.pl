#!/usr/bin/perl 
use strict;
use MakeDB;
use Statistics::Descriptive;


my $DBref = MakeDB();
my @series =();
my $fname = $ARGV[0];
my $begin = $fname;
if ($fname=~/([^\/]+)$/){
    $begin = $1;
    print $begin . "\n";
}

my $dir = $begin . "/";
unless(-e $dir){
    system "mkdir $dir";
}


print STDERR "Parsing divisions for series in $begin\n";
open(IN, $fname) || die "couldnt open $fname\n";
foreach(<IN>){
    chomp;
    my($series, @etc)=split(/\t/);
    if(exists($DBref->{$series})){
	push(@series, $series);
	print STDERR "USING $series\n";
    }else{
	print STDERR "NOT USING $series\n";
    }
}
close IN;

my $series;
my (%x,%y,%z,%d,%gweight,%count);
foreach $series (@series){
    my $annots = $DBref->{$series}{"annots"} . "/dats/SCD$series" . ".csv";
    foreach(`cat $annots`){
	chomp;
	my ($ct, $cell, $time, $none, $global, $local, $blot, $cross, $z, $x, $y, $d, $gweight) = split(/,/);
	next if ($cell eq "cell");
	next if ($x == 0);
	$x{$cell}{$time} += $x;
	$y{$cell}{$time} += $y;
	$z{$cell}{$time} += $z;
	$d{$cell}{$time} += $d;
	$gweight{$cell}{$time} += $gweight;
	$count{$cell}{$time} ++;
    }
}

my $cell;
foreach $cell (sort keys %x){
    my $time;
    foreach $time(sort keys %{$x{$cell}}){
	my $ct =  $count{$cell}{$time};
	if($ct>0){
	    print join(",",$cell . ":" . $time ,$cell, $time, 0,$ct, int($gweight{$cell}{$time}/$ct), 0,0,$z{$cell}{$time}/$ct,int($x{$cell}{$time}/$ct),int($y{$cell}{$time}/$ct),int($d{$cell}{$time}/$ct),int($gweight{$cell}{$time}/$ct)) . "\n";
	}else{
	    print STDERR "CT is $ct  for $cell $time - skipping\n";
	}
    
    }
}
