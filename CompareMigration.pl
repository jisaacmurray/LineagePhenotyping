#!/usr/bin/perl 
use strict;
use Statistics::Descriptive;
use MakeDB;
use Math::Trig;

my $pi = 4*atan2(1,1);
#print STDERR "PI is $pi\n";

my $DBref = MakeDB();


my %series;
my %cells;
my (%x, %y, %z);
my (%r, %theta);
my (%xavg,%yavg,%zavg);


my $filename = $ARGV[0];
open(INPUT, $filename) || die("couldn't open $filename\n");
foreach (<INPUT>){
    chomp;
    my ($series, $endtime, $endcells) = split(/\t/);
#    print STDERR "PROCESSING $series\n";

    my $SCDfile = $DBref->{$series}{"annots"} . "/dats/SCD$series.csv";
    my $AUXINFO = $DBref->{$series}{"annots"} . "/dats/$series" . "AuxInfo.csv";
#    print STDERR "AUXINFO FILE $AUXINFO\n";
    my $xmlfile = $DBref->{$series}{"annots"} . "/dats/$series" . ".xml";

    my $xyres = 0.09;
    my $zres = 0.25;
    my $planeend = 67;
    foreach(`cat $xmlfile`){
	if(/xyRes=\"(.+)\" zRes=\"(.+)\" planeEnd=\"(\d+)/){
#            xyRes=".087  " zRes=".504  " planeEnd="67"/
	    $xyres = $1;
	    $zres = $2;
	    $planeend = $3;
	    $zres = 0.25;
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
#	print "STEPPING THROUGH AUXINFO\n";
	next if(/name/);
	my($name,$slope,$intercept,$xc,$yc,$maj,$min,$ang,$zc,$zslope,$time,$zpixres,$axis) = split(/,/);
#	name,slope,intercept,xc,yc,maj,min,ang,zc,zslope,time,zpixres,axis
#	print;
# 	$xoffset = $xc;
# 	$yoffset = $yc;
# 	$zoffset = $zc;
# 	if($axis =~ /([AaPp])([DdVv])([LlRr])/){
# 	    if($1 eq "P"|| $1 eq "p"){
# 		$xcor = -1;
# 	    }
# 	    if($2 eq "V" || $2 eq "v"){
# 		$ycor = -1;
# 	    }
# 	    if($3 eq "R" || $3 eq "r"){
# 		$zcor = -1;
# 	    }
# 	}
#	
#	$angle = -1 * ($ang/360)*2*$pi;
#	print stderr $ang;
    }

#    print STDERR " $series ANGLE $angle AXIS $xcor $ycor $zcor\n";

    $series{$series} = 1;
    my $xsum;
    my $ysum;
    my $zsum;
    my $avgcount =0;
    my $xmax;
    my $xmin;
    my $zmax;
    my $zmin;
    my $ymax;
    my $ymin;
    open(FILE, $SCDfile) || next;
    foreach(<FILE>){
	chomp;
	my($merge, $cell, $t, $none,$global,$loc,$blot,$cross, $z, $x, $y, $d, $gweight) = split(/,/);
	if($x==0){
	    next;
	}
#	print $_ . "\n";
#	print "$series $cell $t $x $y $z\n";
	$cells{$cell}{$t} = 1;
#	    $x{$cell}{$series} = ($x * $xcor + $xoffset) * $xyres;
#	    $y{$cell}{$series} = ($y * $ycor + $yoffset) * $xyres;

	$x = ($x)*$xcor*$xyres;
	$y = ($y)*$ycor*$xyres;
	$z = ($z)*$zcor*$zres;
	
	

#	my $xnew = (cos($angle)*$x - sin($angle) * $y);
#	my $ynew = (cos($angle)*$y + sin($angle) * $x);
#       don't rotate since SCD file is pre-rotated
	my $xnew = $x;
	my $ynew = $y;

	$xsum += $x{$cell}{$t}{$series} = $xnew;
	$ysum += $y{$cell}{$t}{$series} = $ynew;
	$zsum += $z{$cell}{$t}{$series} = $z;


#	    print "$series $cell x $x y $y z $z     -$angle-     x $xnew y $ynew z $z\n";
	$avgcount++;
	
    }
    if($avgcount>0){
	$xavg{$series}{"all"}=$xsum/$avgcount;
	$yavg{$series}{"all"}=$ysum/$avgcount;
	$zavg{$series}{"all"}=$zsum/$avgcount;
#	print "XAVERAGE FOR $series = $xavg{$series} \n";
#	print "YAVERAGE FOR $series = $yavg{$series} \n";
#	print "ZAVERAGE FOR $series = $zavg{$series} \n";
    }else{
	print STDERR "NO DATA FOR $series\n";
    }
   
    close FILE;
    
}
close INPUT;

if($filename=~/([^\/]+)$/){
    $filename = $1;
    
}

OutputPositions($filename);
exit(0);
# AlignXY(\%z);
# AlignXY(\%x);
# AlignXY(\%y);





sub OutputPositions{
    my $fname = shift;
    my $avgfname = $fname . ".average";
    my %sumsX;
    my %sumsY;
    my %sumsZ;
    my %counts;

    my $dir = $fname . "/";
    
    unless(-e $dir){
	system "mkdir $dir";
    }

    open(OUT, ">$dir/$fname" . "_positions.txt") || die "Couldn't open file $dir/$fname for output\n";
    open(RADIAL, ">$dir/$fname" . "_radialpositions.txt") || die "Couldn't open file $dir/$fname for radial output\n";
    print OUT join("\t",qw(cell time));
    print RADIAL join("\t",qw(cell time));
    my $series;
    foreach $series (sort keys %series){
	print OUT "\t" .join("\t", "X_$series", "Y_$series", "Z_$series");
	print RADIAL "\t" .join("\t", "X_$series", "R_$series", "THETA_$series");
    }
    print OUT "\n";
#   print STDERR "OUTPUTTING Positions $fname\n";
    my $cell;
    foreach $cell (sort keys %cells){
	my $time;
	foreach $time (sort {$a<=>$b} keys %{$x{$cell}}){
	    print OUT join("\t", $cell, $time);
	    print RADIAL join("\t", $cell, $time);

	    my $series;
	    foreach $series (sort keys %series){
		if(exists($x{$cell}{$time}{$series})){
		    my $xnew = SigFigs($x{$cell}{$time}{$series} - $xavg{$series}{"all"},1);
		    my $ynew = SigFigs($y{$cell}{$time}{$series} - $yavg{$series}{"all"},1);
		    my $znew = SigFigs($z{$cell}{$time}{$series} - $zavg{$series}{"all"},1);

		    my $r=sqrt($ynew**2 + $znew**2);
		    my $theta;
		    if($ynew>0){
			$theta = atan($znew/$ynew) ;
		    }elsif($ynew<0){
			if($znew>=0){
			    $theta = atan($znew/$ynew) + $pi;
			}else{
			    $theta = atan($znew/$ynew) - $pi;
			}
		    }else{
			if($znew>0){
			    $theta = $pi/2;
			}elsif($znew<0){
			    $theta = -$pi/2;
			}else{
			    $theta=0;
			}
		    }
#	print STDERR "$cell $t XYZ $xnew $ynew $z XRtheta $xnew $r{$cell}{$t} $theta\n";
		    
		    print OUT "\t" . join("\t", $xnew,$ynew,$znew);
		    print RADIAL "\t" . join("\t", $xnew, SigFigs($r,1),SigFigs($theta,3));
		}else{
		    print OUT "\t" x 3;
		    print RADIAL "\t" x 3;
		}
	    }
	    print OUT "\n";
	    print RADIAL "\n";
	}
    }
    close OUT;
    close RADIAL;
    open (IN, $dir . "/" . $fname . "_positions.txt") || die "Couldn't read $dir/$fname" . "_positions\n";
    open (DEV, ">$dir/$fname" . "_deviations.txt") || die "Couldn't write $dir/$fname deviations\n";
    open (DIST, ">$dir/$fname" . "_distance.txt") || die "Couldn't write $dir/$fname distance\n";
    foreach(<IN>){
	####THIS IS STILL A WORK IN PROGRESS
#	next;
	chomp;
	if(/^cell/){
	    my ($cell, $time, @data) = split(/\t/);
#	    print STDERR "$cell\n";
	    print DEV join("\t", $cell, $time, "avgX", "avgY", "avgZ", @data) . "\n";

	    print DIST join("\t", $cell, $time);
	    for(my $i=0;$i<@data;$i+=3){
		print DIST "\t";		
		if($data[$i]=~/[XYZ](_.+)/){
		    print DIST "DIST_$1";
		}
	    }
	    print DIST "\n";

	}else{
	    my ($cell,$time, @data) = split(/\t/);
#	    print STDERR "$cell $time\n";
	    my $count=0;
	    my $ySum=0;
	    my $zSum=0;
	    my $xSum=0;
	    
	    for (my $i=0;$i<@data;$i+=3){
		my $x = $data[$i];
		my $y = $data[$i+1];
		my $z = $data[$i+2];
		if($x=~/\d/){
		    $count++;
		    $xSum += $x;
		    $ySum += $y;
		    $zSum += $z;
		}
	    }

	    if($count>0){
#		print "$cell $time $xSum $ySum $zSum $count\n";
		print DEV join("\t", $cell, $time);
		print DIST join("\t", $cell, $time);
		my $xavg = $xSum/$count;
		my $yavg = $ySum/$count;
		my $zavg = $zSum/$count;

		print DEV "\t" . join("\t", SigFigs($xavg,1), SigFigs($yavg,1), SigFigs($zavg,1));
		
		for (my $i=0;$i<@data;$i+=3){
		    my $x = $data[$i];
		    my $y = $data[$i+1];
		    my $z = $data[$i+2];
		    if($x=~/\d/){
			print DEV "\t" . join("\t", SigFigs($x-$xavg,1),
					  SigFigs($y-$yavg,1),
					  SigFigs($z-$zavg, 1));
			print DIST "\t" . sqrt(($x-$xavg)**2 + ($y-$yavg)**2 + ($z-$zavg)**2);
		    }else{
			print DEV "\t" x 3;
			print DIST "\t";
		    }
		}
		print DEV "\n";
		print DIST "\n";
	    }else{
#		print STDERR "SKIPPING $cell $time - COUNT $count <= 0\n";
	    }
	}
    }
    close DEV;
    close DIST;

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
	    if ($hash->{$cell}{$series} =~/\d/){
		$stat->add_data($hash->{$cell}{$series});
	    }
	    if($hash->{$cell}{$series} =~/\d/ &&
	       $hash->{$cell}{$reference} =~/\d/){
		$stat1->add_data($hash->{$cell}{$series});
		$stat2->add_data($hash->{$cell}{$reference});
	    }
	}
	
	my $mean = $stat->mean();
	my ($q, $r, $m, $rms) = $stat1->least_squares_fit($stat2->get_data());
	next if ($m == 0);
 	foreach $cell (sort keys %cells){
	    if ($hash->{$cell}{$series} =~/\d/){
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
    open(OUT, ">$dir/$fname");
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
    open(OUT, ">$dir/$fname");
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
sub SigFigs{
    my ($value, $figs) = @_;
    return int($value * (10**$figs))/(10**$figs);
}
