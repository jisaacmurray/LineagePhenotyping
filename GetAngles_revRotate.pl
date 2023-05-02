#!/usr/bin/perl -w
use strict;
use Statistics::Descriptive;
use MakeDB;

my $pi = 4*atan2(1,1);
#print STDERR "PI is $pi\n";


my %refX;
my %refY;
my %refZ;
#my $referenceAngles = "/Users/jmurr/Dropbox/LR_TF_Function/ceh-36/Over500WT/Over500WT.angles_aligned.txt.average";
my $referenceAngles = "Richard_et_al_plus_comma_WT/Richard_et_al_plus_comma_WT.angles_aligned.txt.average";
open(REFANG, $referenceAngles) || die "Couldn't open $referenceAngles\n";
foreach(<REFANG>){
    chomp;
    my($cell, $x, $y, $z) = split(/\t/);
    if($x=~/\d/){
	$refX{$cell} = $x;
	$refY{$cell} = $y;
	$refZ{$cell} = $z;
    }
}

close REFANG;

my $DBref = MakeDB();

my %series;
my %cells;
my (%x, %y, %z);

my $filename = $ARGV[0];
open(FILE, $filename) || die("couldn't open $filename\n");
foreach (<FILE>){
    chomp;
    my ($series, $endtime, $endcells) = split(/\t/);
#    print STDERR "PROCESSING $series\n";

    my $CDfile = $DBref->{$series}{"annots"} . "/dats/CD$series.csv";
    my $AUXINFO = $DBref->{$series}{"annots"} . "/dats/$series" . "AuxInfo.csv";
#    print STDERR "AUXINFO FILE $AUXINFO\n";
    my $xmlfile = $DBref->{$series}{"annots"} . "/dats/$series" . ".xml";
    open(FILE, $CDfile) || next;

    my $xyres = 0.09;
    my $zres = 0.5;
    my $planeend = 67;
    foreach(`cat $xmlfile`){
	if(/xyRes=\"(.+)\" zRes=\"(.+)\" planeEnd=\"(\d+)/){
#            xyRes=".087  " zRes=".504  " planeEnd="67"/
	    $xyres = $1;
	    $zres = $2;
	    $planeend = $3;
	}else{
#	    print "NOTRESLINE: $_";
	}
    }
#    print STDERR "$series XYRES $xyres ZRES $zres\n";


    my $xcor=1;
    my $ycor=1;
    my $zcor=1;
    my $xoffset = 0;
    my $yoffset = 0;
    my $zoffset = 0;
    my $angle = 0;
    foreach(`cat $AUXINFO`){
	chomp;
#	print "STEPPING THROUGH AUXINFO\n";
	next if(/name/);
	my($name,$slope,$intercept,$xc,$yc,$maj,$min,$ang,$zc,$zslope,$time,$zpixres,$axis) = split(/,/);
#	name,slope,intercept,xc,yc,maj,min,ang,zc,zslope,time,zpixres,axis
#	print;
	print STDERR "$series $axis\n";
	if($axis =~ /([AaPp])([DdVv])([LlRr])/){
	    if($1 eq "P"|| $1 eq "p"){
		$xcor = -1;
		$xoffset = 712;
	    }
	    if($2 eq "V" || $2 eq "v"){
		$ycor = -1;
		$yoffset = 512;
	    }
	    if($3 eq "R" || $3 eq "r"){
		$zcor = -1;
		$zoffset = $planeend;
	    }
	}
	
	$angle =  ($ang/360)*2*$pi;
	print STDERR $ang;
    }

    print STDERR " $series ANGLE $angle AXIS $xcor $ycor $zcor\n";

    $series{$series} = 1;
    my $lastCell = 0;
    foreach(<FILE>){
	chomp;
	next if(/^cell/);
	my($merge, $cell, $t, $none,$global,$loc,$blot,$cross, $z, $x, $y, $d, $gweight) = split(/,/);

	$cells{$cell} = 1;
	
	if($lastCell ne $cell){
#	    $x{$cell}{$series} = ($x * $xcor + $xoffset) * $xyres;
#	    $y{$cell}{$series} = ($y * $ycor + $yoffset) * $xyres;
	    
	    my $xnew = cos($angle)*$x - sin($angle) * $y;
	    my $ynew = cos($angle)*$y + sin($angle) * $x;
	    
	    $x{$cell}{$series} = ($xnew * $xcor + $xoffset) * $xyres;
	    $y{$cell}{$series} = ($ynew * $ycor + $yoffset) * $xyres;

	    $z{$cell}{$series} = ($z * $zcor + $zoffset) * $zres;
#	    print "$series $cell x $x y $y z $z     -$angle-     x $xnew y $ynew z $z\n";

	}
	$lastCell = $cell;
    }
   
    close FILE;
    
}

if($filename=~/([^\/]+)$/){
    $filename = $1;
}


OutputAngles($filename . ".angles_unaligned.txt", $filename);

AlignXY(\%z);
AlignXY(\%x);
AlignXY(\%y);

OutputAngles($filename . ".angles_aligned.txt", $filename);



sub OutputAngles{
    my $fname = shift;
    my $dir = shift;
    my $avgfname = $fname . ".average";


    my %sumsX;
    my %sumsY;
    my %sumsZ;
    my %counts;
    $dir = $dir . "/";
    unless(-e $dir){
	system "mkdir $dir\n";
    }
    open(OUT, ">$dir/$fname") || die "Couldn't open file $dir/$fname for output\n";
    print OUT join("\t",qw(series parent daughter1 x1 y1 z1 daughter2 x2 y2 z2 dx dy dz refx refy refz dot)) . "\n";
    print STDERR "OUTPUTTING ANGLES $fname\n";
    my $series;
    foreach $series (sort keys %series){
#	print STDERR "SERIES $series\n";
	my $cell;

	foreach $cell (sort keys %cells){

	    my @daughters;


	    if($cell eq "EMS"){
		@daughters = ("E", "MS");
	    } elsif ($cell eq "P0"){
		@daughters = ("AB", "P1");
	    } elsif ($cell eq "P1"){
		@daughters = ("EMS", "P2");
	    } elsif ($cell eq "P2"){
		@daughters = ("C", "P3");
	    } elsif ($cell eq "P3"){
		@daughters = ("D", "P4");
	    } elsif ($cell eq "P4"){
		@daughters = ("Z2", "Z3");
	    }

	    foreach(sort keys %cells){
		if(/$cell.$/){
#		    print STDERR "$cell daughters include $_\n";	
		    push(@daughters, $_);
		}else{
#		    print STDERR "$_ not a daughter of $cell\n";
		}
	    }
	    if(defined($daughters[0]) && defined($daughters[1]) && defined($x{$daughters[0]}{$series}) && defined($x{$daughters[1]}{$series}) && 
	       @daughters >0 && $x{$daughters[0]}{$series} =~/\d/ && $x{$daughters[1]}{$series} =~/\d/){
		
		print OUT join("\t", $series, $cell);
		foreach(@daughters){
		    print OUT "\t" . join("\t", $_, $x{$_}{$series},$y{$_}{$series},$z{$_}{$series});
		}
		my $dx = $x{$daughters[1]}{$series} - $x{$daughters[0]}{$series};
		my $dy = $y{$daughters[1]}{$series} - $y{$daughters[0]}{$series};
		my $dz = $z{$daughters[1]}{$series} - $z{$daughters[0]}{$series};
		my $mag= sqrt($dx**2 + $dy**2 + $dz**2);
		unless($mag==0){
		    $dx /= $mag;
		    $dy /= $mag;
		    $dz /= $mag;
		    $sumsX{$cell} += $dx;
		    $sumsY{$cell} += $dy;
		    $sumsZ{$cell} += $dz;
		    $counts{$cell} ++;
		}

		

		print OUT "\t" . join("\t", $dx, $dy, $dz);

		if(!defined($refX{$cell})){
		    $refX{$cell} = 1;
		    $refY{$cell} = 0;
		    $refZ{$cell} = 0;
		}

		print OUT "\t" . join("\t", $refX{$cell}, $refY{$cell}, $refZ{$cell}, $dx * $refX{$cell} + $dy*$refY{$cell} + $dz*$refZ{$cell});
		
		print OUT "\n";
	    }
	}
	
    }
    close OUT;
    open(AVG, ">$dir/$avgfname") || die "Couldn't open file $dir/$avgfname for output\n";
    foreach(sort keys %counts){
	print AVG join("\t", $_, 
		   $sumsX{$_}/$counts{$_},
		   $sumsY{$_}/$counts{$_},
		   $sumsZ{$_}/$counts{$_}) . "\n";
    }
    close AVG;
}


sub Zcorrect{
    my($hash1, $z, $factor, $center) = @_;
    my $cell;
    my $series;
    foreach $cell (sort keys %cells){
	foreach $series (sort keys %series){
	    if($hash1->{$cell}{$series} =~/\d/ &&
	       $z->{$cell}{$series} =~/\d/){
		$hash1->{$cell}{$series} = $hash1->{$cell}{$series} / (1-$factor*($z->{$cell}{$series} - $center));

	    }
	    else{
		$hash1->{$cell}{$series} = "";
	    }
	}
    }
}


sub AlignXY{
    my $hash = shift;
    my $reference = "20110316_mir-1_L2";

    my ($series, $cell);


    foreach $series (sort keys %series){
	my $stat = Statistics::Descriptive::Full->new();
	my $stat1 = Statistics::Descriptive::Full->new();
	my $stat2 = Statistics::Descriptive::Full->new();
	
	foreach $cell (sort keys %cells){
	    next if(!defined($hash->{$cell}{$series}));
	    if ($hash->{$cell}{$series} =~/\d/){
		$stat->add_data($hash->{$cell}{$series});
	    }
	    if(defined($hash->{$cell}{$reference}) &&
	       $hash->{$cell}{$series} =~/\d/ &&
	       $hash->{$cell}{$reference} =~/\d/){
		$stat1->add_data($hash->{$cell}{$series});
		$stat2->add_data($hash->{$cell}{$reference});
	    }
	}
	
	my $mean = $stat->mean();
	my ($q, $r, $m, $rms) = $stat1->least_squares_fit($stat2->get_data());
	next if ($m == 0);
 	foreach $cell (sort keys %cells){
	    if (defined($hash->{$cell}{$series}) && $hash->{$cell}{$series} =~/\d/){
		$hash->{$cell}{$series} = ($hash->{$cell}{$series} - $mean)/$m;
		
	    }
	}
    }
  

}





sub Normalize{
    my $hash = shift;

    my ($series, $cell);


    foreach $series (sort keys %series){
	my $stat = Statistics::Descriptive::Full->new();

	foreach $cell (sort keys %cells){
	    if ($hash->{$cell}{$series} =~/\d/){
		$stat->add_data($hash->{$cell}{$series});
	    }
	}
	
	my $mean = $stat->mean();

 	foreach $cell (sort keys %cells){
	    if ($hash->{$cell}{$series} =~/\d/){
		if($mean == 0){
		    $hash->{$cell}{$series} = "";
		}else{
		    $hash->{$cell}{$series} = ($hash->{$cell}{$series} / $mean);
		}
		
	    }
	}
    }
}


sub TopAlign{
    my $hash = shift;

    my ($series, $cell);


    foreach $series (sort keys %series){
	my $stat = Statistics::Descriptive::Full->new();


	foreach $cell (sort keys %cells){
	    if ($hash->{$cell}{$series} =~/\d/){
		$stat->add_data($hash->{$cell}{$series});
	    }
	}
#	print STDERR "Calculating 10th percentile for $cell $series\n";
	my $tenth = $stat->percentile(10);

 	foreach $cell (sort keys %cells){
	    if ($hash->{$cell}{$series} =~/\d/){
		$hash->{$cell}{$series} = ($hash->{$cell}{$series} - $tenth);
		
	    }
	}
    }
    
}



sub Correlate{
    my ($fname, $hash1, $hash2) = @_;
    open(OUT, ">$fname");
    print OUT "Cell\tComparison\tmean1\tstdev1\tn1\tmean2\tstdev2\tn2\tq\tm\tr\trms\n ";
    print STDERR "Ouputting Cor to $fname\n";
    # Data
    my ($series, $cell);
    foreach $cell (sort keys %cells){
	my (@ar1, @ar2);
	
	foreach $series (sort keys %series){
	    if($hash1->{$cell}{$series} =~/\d/ &&
	       $hash2->{$cell}{$series} =~/\d/){
		push(@ar1,  $hash1->{$cell}{$series});
		push(@ar2,  $hash2->{$cell}{$series});
	    }
	}
	
	

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@ar1);
	my ($q, $m, $r, $rms) = $stat->least_squares_fit(@ar2);
	my $stat2 = Statistics::Descriptive::Full->new();
	$stat2->add_data(@ar2);

	print OUT join ("\t",$cell, $fname, 
			$stat->mean(), $stat->standard_deviation(), $stat->count(),
			$stat2->mean(), $stat2->standard_deviation(), $stat2->count(),
			$q, $m, $r, $rms) . "\n";

    }

    close OUT;
}




sub Output{

    my ($fname, $hash) = @_;
    print STDERR "Outputting $fname\n";
    open(OUT, ">$fname");
    my ($series, $cell);

    #print headers
    print OUT "Cell";
    foreach $series (sort keys %series){
	print OUT "\t$series";
    }
    print "\n";
    
    #print Data
    foreach $cell (sort keys %cells){
	print OUT $cell;
	foreach $series (sort keys %series){
	    print OUT "\t" . $hash->{$cell}{$series};
	}
	print OUT "\n";
    }
    close OUT;
}
