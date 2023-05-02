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


my %min;
my %max;
my @series = ($sulston);

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


open(OUT, ">$dir/$begin" . "CellsVsTime.tsv") || die "Couldn't write to $dir/$begin" . "CellsVsTime.tsv\n";
print OUT join("\t", qw(Series Minutes TPs nCells AB E MS C D)) . "\n";

my $series;
foreach $series (@series){



#    print STDERR "READING ANNOTS FOR $series\n";
    my $annots = $DBref->{$series}{"annots"} . "/dats/CD$series" . ".csv";
    my $timefile = $DBref->{$series}{"annots"} . "/dats/TIME$series" . ".csv";
    
    my (%AB, %E, %MS, %C, %D, %total);

    if(-e $annots && (-e $timefile || $series eq $sulston)){
	my %timelookup;
	if(-e $timefile){
	    foreach(`cat $timefile`){
		my ($seriesDep, $tp, $seconds, $delta) = split(/\t/);
		$tp = int($tp);
		$timelookup{1*$tp} = (int($seconds/6))/10;
#		print "$tp $timelookup{$tp}\n";
	    }
	}elsif($series ne $sulston){
	    next;
	    #Skip data acquisition unless 
	}
	print STDERR "READING ANNOTS $annots FOR SERIES $series\n";
	foreach(`cat $annots`){
#	    print;
	    chomp;
	    my ($ct, $cell, $time, $none, $global, $local, $blot, $cross, $z, $x, $y, $d, $gweight) = split(/,/);
	    next if ($cell eq "cell");

#	    print "$cell $time\n";
	    if($series ne $sulston){
		if(defined($timelookup{$time})){
		    $time = $timelookup{$time};
		}
#		print "$cell $time\n";
	    }


	    #UPDATE COUNTS OF CELLS PER EACH LINEAGE vs TIME
	    if($cell=~/^AB/){
		$AB{$time}++;
	    }elsif($cell=~/^E/ && $cell ne "EMS"){
		$E{$time}++;
	    }elsif($cell=~/^MS/){
		$MS{$time}++;
	    }elsif($cell=~/^C/){
		$C{$time}++;
	    }elsif($cell=~/^D/){
		$D{$time}++;
	    }
	    $total{$time}++;
	    
	    


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
	    

	    #STORE MIN AND MAX (OBSERVED) TIME FOR EACH CELL IN HASH
	    if(!exists($min{$cell}{$series}) || $min{$cell}{$series} > $time){
		$min{$cell}{$series} = $time;
#		print "MIN $cell $series $time\n";
	    }
	    if(!exists($max{$cell}{$series}) || $max{$cell}{$series} < $time){
		$max{$cell}{$series} = $time;
#		print "MAX $cell $series $time\n";
	    }
	}
    }else{
	print STDERR "Couldn't find $annots - SKIPPING $series\n";
    }
    
    my $tp=1;
    #OUTPUT CELLS VS TIME TO FILE
    foreach(sort{$a<=>$b} keys %total){
	next unless ($_=~/\d/);
	print OUT join("\t", $series, $_, $tp++, map {$_=~/\d/?$_:0} ($total{$_},$AB{$_},$E{$_},$MS{$_},$C{$_},$D{$_})) . "\n";
    }
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
open(CC, (">$dir/$begin" . "CCLength.tsv")) || die "couldnt open CC\n";;
open(DivTime, (">$dir/$begin" . "DivTime.tsv"))|| die "couldnt open DIVTIME\n";;

print CC "Cell\t" . join("\t", @series) . "\n";
print DivTime "Cell\t" . join("\t",@series) . "\n";
my %ratios;
#cycle through all cells, calculate rates for internal branches and print raw cell cycle and cell division times
foreach $cell(sort keys %min){ #All cells identified in ANY series
#    print STDERR "PARSING DIVTIMES FOR $cell\n";
    print CC $cell;
    print DivTime $cell;
#    print STDERR $cell;
    
    foreach $series(@series){
	my $cc =  ($max{$cell}{$series} - $min{$cell}{$series});
	my $sulstoncc = ($max{$cell}{$sulston} - $min{$cell}{$sulston});
	
	#PRINT CCTIME OR MINCCTIME DEPENDING ON WHETHER THIS CELL DIVIDED AND BIRTH WAS OBSERVED
	#ALSO STORE RATIO IF THIS IS A GOOD DIVISION IN BOTH SULSTON AND THIS SERIES
	if(!HasProgeny($cell, $series, \%progeny) || !HasParent($cell,$series,\%parents)){
#	    print CC "\t>" . SigFigs($cc,1);	    
	    print CC "\tNA";
	}else{
	    print CC "\t" . SigFigs($cc,1);
	    if($sulstoncc>0){
		$ratios{$series}{$cell} = $cc/$sulstoncc;
	    }

	}
	#PRINT DIVTIME OR MINDIVTIME DEPENDING ON WHETHER THIS CELL DIVIDED
	if(!HasProgeny($cell, $series, \%progeny)){
#	    print DivTime "\t>" .  SigFigs($max{$cell}{$series},1);
	    print DivTime "\tNA";
	}else{
	    print DivTime "\t" .  SigFigs($max{$cell}{$series},1);
	}
	#missing cells/divisions from series, very late divisions in sulston, etc
    }
    print CC"\n";
    print DivTime"\n";
}

close CC;
close DivTime;



#calculate division rate slope and offset relative to offset

my %slopes;
my %offsets;
my %lineageSlopes;
foreach $series(@series){
    my $count = 0;
    my $sum = 0;
    my $ratioStat = Statistics::Descriptive::Full->new();

    my %lineageRateBins;
    my %sulstonRateBins;
    my @seriesTimes;
    my @sulstonTimes;
    foreach $cell(sort keys %min){
	#KEEP CELLS FOR ANALYSIS IF 1)BIRTH WAS OBSERVED 2)DIVISION WAS OBSERVED 3)DIVIDES IN SULSTON
	if(HasParent($cell, $series, \%parents) && HasProgeny($cell,$series,\%progeny) && HasProgeny($cell, $sulston, \%progeny)){

	    #update 4 metrics for series slope: 1)average ratio (sum, count), 2)average ratio (ratioStat),
	    #                                   3)divtime metrics 4)lineage-specific rate bins
	    $sum+= $ratios{$series}{$cell};
	    $count++;
	    $ratioStat->add_data($ratios{$series}{$cell});
	    
	    push(@seriesTimes, $max{$cell}{$series});
	    push(@sulstonTimes, $max{$cell}{$sulston});
	    
	    foreach(@founders){
		if($cell =~/$_[apdvlr]*/){
		    push(@{$lineageRateBins{$_}}, $max{$cell}{$series});
		    push(@{$sulstonRateBins{$_}}, $max{$cell}{$sulston});
		}
	    }
	}
    }

    next unless ($count>0);

    #THIS COMMENTED OUT CODE CALCULATES SLOPE BY ONE OF THE TWO CC-based METHODS   
    #    $slopes{$series} = $ratioStat->mean();
    #    $slopes{$series} = $sum/$count;

    #NEW STAT: BASED ON DIVISION TIMES
    my $divtimeStat = Statistics::Descriptive::Full->new();
    $divtimeStat->add_data(@seriesTimes);
    my ($q,$m,$r,$rms) = $divtimeStat->least_squares_fit(@sulstonTimes);
    $slopes{$series} = $m;
    if($slopes{$series} ==0){
	print STDERR "$series HAS ZERO SLOPE - WILL FAIL!\n";
    }
    $offsets{$series} = $q;
    

    foreach(@founders){
	my $lineageStat =  Statistics::Descriptive::Full->new();
	$lineageStat->add_data(@{$lineageRateBins{$_}});
	my ($qLin,$mLin,$rLin,$rmsLin) = $lineageStat->least_squares_fit(@{$sulstonRateBins{$_}});
	if($mLin == 0){
	    $lineageSlopes{$series}{$_} = $slopes{$series};
	    print STDERR "FAILED TO FIND SLOPE FOR $series $_ - USING AVERAGE\n";
	}else{
	    $lineageSlopes{$series}{$_} = $mLin;
	}
    }

    print STDERR join(" ", "$series CC RATIO MEANS: ", $ratioStat->mean(), $sum/$count, "DIVTIME REGRESSION SLOPE/INTERCEPT:", $m, $q) . "\n";


    #FORMERLY:calculate offset by correcting times and scaling relative to fixed landmark division ABalaa division at 103 minutes (Sulston)
    # and use CC ratio with this offset to correct series
}

open(CCNORM, (">$dir/$begin" . "CCLengthNorm.tsv"));
open(DivTimeNORM, (">$dir/$begin" . "DivTimeNorm.tsv"));
open(CCLINNORM, (">$dir/$begin" . "CCLinNorm.tsv"));
open(CCNORMTERMINAL, (">$dir/$begin" . "CCLengthMinTerminal.tsv"));
#print CCNORM "Cell\t" . join("\t",@series, qw(count mean stdev CV)) . "\n";
#print DivTimeNORM "Cell\t" . join("\t", @series, qw(count mean stdev CV)) . "\n";
#print CCLINNORM "Cell\t" . join("\t",  @series, qw(count mean stdev CV)) . "\n";
print CCNORM "Cell\t" . join("\t",@series) . "\n";
print DivTimeNORM "Cell\t" . join("\t", @series) . "\n";
print CCLINNORM "Cell\t" . join("\t",  @series) . "\n";
print CCNORMTERMINAL "Cell\t" . join("\t",  @series) . "\n";

open(STATS,  (">$dir/$begin" . "STATS.tsv"));
print STATS join("\t", qw(series offset slope), @founders) . "\n";
foreach $series (@series){
    print STATS join("\t", $series, SigFigs($offsets{$series},1), SigFigs($slopes{$series},3), map({SigFigs($lineageSlopes{$series}{$_},3)} (@founders)) ) . "\n";
}

close STATS;

foreach $cell(sort keys %min){
    next unless ($cell=~/\w+/);
    print CCNORM $cell;
    print CCLINNORM $cell;
    print DivTimeNORM $cell;
    print CCNORMTERMINAL $cell;
    my $CCnormstat =  Statistics::Descriptive::Full->new();
    my $CClinnormstat =  Statistics::Descriptive::Full->new();
    my $DivTimeNormStat =  Statistics::Descriptive::Full->new();

    
    foreach $series(@series){
	
	my $cc = 0;
	if($slopes{$series} != 0){
	    $cc = ($max{$cell}{$series} - $min{$cell}{$series})/$slopes{$series};
#	    print STDERR "CC $cc $series $cell $max{$cell}{$series} - $min{$cell}{$series}\n";
	}
	my $sulstoncc = ($max{$cell}{$sulston} - $min{$cell}{$sulston});
	if(!HasParent($cell,$series, \%parents) || !HasProgeny($cell,$series, \%progeny)){
#	    print CCNORM "\t>" . SigFigs($cc,1);
	    print CCNORM "\tNA";
	    
	    if($cc>0){
		print CCNORMTERMINAL "\t" . SigFigs($cc,1);
	    }else{
		print CCNORMTERMINAL "\tNA";
	    }
	    #LinNorm leaves out missing values still
	    print CCLINNORM "\tNA";

	}else{
	    print CCNORM "\t" . SigFigs($cc,1);
	    unless($series eq $sulston){
		$CCnormstat->add_data($cc);
	    }
	    print CCNORMTERMINAL "\tNA";
	    my $printed = 0;
	    foreach(@founders){
		if($cell=~/$_[apdvlr]*/){
		    $printed =1;
#		    print STDERR "USING SLOPE $lineageSlopes{$series}{$_} for $cell founder $_\n";
		    unless($lineageSlopes{$series}{$_} >0){
			$lineageSlopes{$series}{$_} = 1;
			print STDERR "LINEAGE SLOPE IS $lineageSlopes{$series}{$_} FOR SERIES $series FOUNDER $_ - RESETTING TO 1\n";
		    }
		    print CCLINNORM "\t" . SigFigs((($max{$cell}{$series} - $min{$cell}{$series}) / $lineageSlopes{$series}{$_}),1);
		    unless($series eq $sulston){
			$CClinnormstat->add_data(($max{$cell}{$series} - $min{$cell}{$series}) / $lineageSlopes{$series}{$_});
		    }
#		    print STDERR "FOUNDER FOR $cell is $_\n";
		}
	    }
	    if($printed==0){
		print CCLINNORM "\t" . SigFigs($cc,1);
#		print STDERR "Couldn't find founder for $cell\n";
	    }

	}
	my $thisslope = $slopes{$series};
	if($thisslope == 0){
	    $thisslope = 1;
	    print STDERR "SLOPE FAILED FOR $series\n";
	}
	if(!HasProgeny($cell,$series, \%progeny)){
#	    print DivTimeNORM "\t>" .  SigFigs(($max{$cell}{$series}-$offsets{$series})/$slopes{$series},1);
	    #BELOW IS GREATER THAN FORMAT
#	    print DivTimeNORM "\t>" .  SigFigs(($max{$cell}{$series}-$offsets{$series})/$thisslope,1);
	    print DivTimeNORM "\tNA";
	}else{
#	    print DivTimeNORM "\t" .  SigFigs(($max{$cell}{$series}-$offsets{$series})/$slopes{$series},1);
	    print DivTimeNORM "\t" .  SigFigs(($max{$cell}{$series}-$offsets{$series})/$thisslope,1);
	    unless($series eq $sulston){
#		$DivTimeNormStat->add_data(($max{$cell}{$series}-$offsets{$series})/$slopes{$series},1);
		$DivTimeNormStat->add_data(($max{$cell}{$series}-$offsets{$series})/$thisslope,1);
	    }
	}
    }
    

    # if($CCnormstat->count()>2){
    # 	print CCNORM "\t" . join("\t", $CCnormstat->count(),$CCnormstat->mean(), $CCnormstat->standard_deviation(), $CCnormstat->standard_deviation()/$CCnormstat->mean()) .  "\n";
    # }else{
    # 	print CCNORM "\t\t\t\n";
    # }


    # if($DivTimeNormStat->count()>2){
    # 	print DivTimeNORM "\t" . join("\t", $DivTimeNormStat->count(),$DivTimeNormStat->mean(), $DivTimeNormStat->standard_deviation(), $DivTimeNormStat->standard_deviation()/$DivTimeNormStat->mean()) .  "\n";
    # }else{
    # 	print DivTimeNORM "\t\t\t\n";
    # }


    # if($CClinnormstat->count()>2){
    # 	print CCLINNORM "\t" . join("\t",$CClinnormstat->count(),  $CClinnormstat->mean(), $CClinnormstat->standard_deviation(), $CClinnormstat->standard_deviation()/$CClinnormstat->mean()) .  "\n";
    # }else{
    # 	print CCLINNORM "\t\t\t\n";
    # }
    print CCNORM "\n";
    print DivTimeNORM "\n";
    print CCLINNORM "\n";
    print CCNORMTERMINAL "\n";
}


close CCNORM;
close CCLINNORM;
close DivTimeNORM;
close CCNORMTERMINAL;



sub HasParent{
    my($cell, $series,$parentref) = @_;
    my $parent = $parentref->{$cell};
    if(exists($min{$parent}{$series})){
	return 1;
    }else{
#	print "No parent $parent for $cell in $series\n";
	return 0;
    }
    
}

sub HasProgeny{
    my($cell, $series,$progref) = @_;
    my $prog = $progref->{$cell};
    if(exists($min{$prog}{$series})){
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
