#!/usr/bin/perl 
use strict;
use MakeDB;
use Statistics::Descriptive;


my $DBref = MakeDB();
my $sulston = "20081128_sulston";
my %parents =(
    "AB" => "P0",
    "P1" => "P0",
    "EMS" => "P1",
    "P2" => "P1",
    "E" => "EMS",
    "MS" => "EMS",
    "C" => "P2",
    "P3" => "P2",
    "D" => "P3",
    "P4" => "P3",
    "Z2" => "P4",
    "Z3" => "P4"
    );

my @founders=qw(ABala ABalp ABara ABarp ABpla ABplp ABpra ABprp E MS C D);


my %blot;
my @series;

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
foreach $series (@series){



#    print STDERR "READING ANNOTS FOR $series\n";
    my $annots = $DBref->{$series}{"annots"} . "/dats/SCA$series" . ".csv";
    
    my (%AB, %E, %MS, %C, %D, %total);

    if(-e $annots){
	my %timelookup;
	print STDERR "READING ANNOTS $annots FOR SERIES $series\n";
	foreach(`cat $annots`){
#	    print;
	    chomp;
	    my ($ct, $cell, $time, $none, $global, $local, $blot, $cross, $z, $x, $y, $d, $gweight) = split(/,/);
	    next if ($cell eq "cell");

	    


	    #UPDATE PARENT MATRIX USING STANDARD NAMING FOR NOVEL CELLS
	    if(!exists($parents{$cell})){
		my $tmp = $cell;
		chop $tmp;
		$parents{$cell} = $tmp;
#		print "PARENT OF $cell IS $tmp\n";
		
	    }
#	    else{
#		print "HARDCODED PARENT OF $cell IS $parents{$cell}\n";
#	    }
	    

	    #STORE BLOT FOR EACH CELL IN HASH
		$blot{$cell}{$series} = $blot;
	    
	}
    }else{
	print STDERR "Couldn't find $annots - SKIPPING $series\n";
    }
    
    my $tp=1;
    #OUTPUT CELLS VS TIME TO FILE
}
close OUT;

my %progeny;
my $cell;
foreach $cell (sort keys %parents){
    #for this bookeeping purpose only need to track one daughter
    # - will use to keep track of whether a cell divides in a given series
    $progeny{$parents{$cell}} = $cell;
}


#OUTPUT RAW CC AND DIVTIME VALUES, STORE RATIOS
open(BLOT, (">$dir/$begin" . "Blot.tsv"))|| die "couldnt open BLOT\n";;

print BLOT "Cell\t" . join("\t",@series) . "\n";
foreach $cell(sort keys %blot){ #All cells identified in ANY series
#    print STDERR "PARSING DIVTIMES FOR $cell\n";
    print BLOT $cell;
#    print STDERR $cell;
    
    foreach $series(@series){
	print BLOT "\t" . $blot{$cell}{$series};
    }
    print BLOT"\n";
}

close BLOT;

sub HasParent{
    my($cell, $series,$parentref) = @_;
    my $parent = $parentref->{$cell};
    if(exists($blot{$parent}{$series})){
	return 1;
    }else{
#	print "No parent $parent for $cell in $series\n";
	return 0;
    }
    
}

sub HasProgeny{
    my($cell, $series,$progref) = @_;
    my $prog = $progref->{$cell};
    if(exists($blot{$prog}{$series})){
#	print "PROGENY FOR $cell in $series!\n";
	return 1;
    }else{
#	print "No progeny $prog for $cell in $series\n";
	return 0;
    }
    
}

sub SigFigs{
    my ($value, $figs) = @_;
    return int($value * (10**$figs))/(10**$figs);
}
