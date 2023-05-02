#!/usr/bin/perl
use strict;
#use MakeDB;
my $seriesList=$ARGV[0];
#my $db=MakeDB();

if(!-e $seriesList){
    die("Couldn't open $seriesList");
}


my @series;
foreach(`cat $seriesList`){
    chomp;
#    print;
    my ($series,@etc)= split(/\t/);
    push(@series, $series);
}


#1) Get average div time in sulston rate-normalized space for each cell
#2) Convert XYZ coordinates of each embryo's cell to common frame of reference using AuxInfo 
#   a)Invert/Flip by ADL/coordinate system to ADL - use auxinfo xc/yc/zc 
#   b)Stretch X,Y by ellipse diameters to match standard ellipse (hard code based on average ellipse)
#   c)Stretch Z by ?  First, try ignoring.  Also consider using top/bottom percentiles?
#IS ZPIXRES THE PROBLEM IN SCD CREATION?

my(%birth, %division,%parent,%daughter);
foreach(`cat SupplementalTable2_DivisionTimes.txt`){
    chomp;
    my($cell,$raw,$norm,$parent, $birth,@etc)=split(/\t/);
#    print "$cell $norm\n";;
    $birth{$cell}=$birth;
    $division{$cell}=$norm;
    $parent{$cell}=$parent;
    #only store one daughter - not important to have both for this application
    $daughter{$parent} = $cell;
}






my %CDs;
my %data;
my $CDheader;
foreach(@series){
#UNCOMMENT THIS AND THE MAKEDB REFS ABOVE TO GRAB FRESH CD FILES
#     $CDs{$_}= $db->{$_}{"annots"} . "/dats/CD" . $_ . ".csv";
#     if(-e  $CDs{$_}){
# 	print "$_ CD FOUND\n";
# 	`cp $CDs{$_} CDs`;
# 	my $AuxInfo=$db->{$_}{"annots"} . "/dats/" . $_ . "AuxInfo.csv";
# 	`cp $AuxInfo  AuxInfos`;
#     }else{
# 	print "$_ CD NOT FOUND\n";
#     }

    my $AuxInfo=$_ . "AuxInfo.csv";
    my($name,$slope,$intercept,$xc,$yc,$maj,$min,$ang,$zc,$zslope,$time,$zpixres,$axis);

    foreach(`cat AuxInfos/$AuxInfo`){
	chomp;
	next if(/name/);
	($name,$slope,$intercept,$xc,$yc,$maj,$min,$ang,$zc,$zslope,$time,$zpixres,$axis) = split(/,/);
    }

    my $CDFile="CD" . $_ . ".csv";
    my $CDOutFile="ACD" . $_ . ".csv";
#    print "$name $xc $yc $zc $axis $maj $min\n";
    open (OUT, ">CDs/$CDOutFile");
    my %loc;
    my @ratios;
    my $prev="";
    my $birthTime=0;
    my $prevTime =0;
    my %thisMaxTime;
    my %thisMinTime;

    my $PI = 4 * atan2(1, 1);
    #rotation angle
    my $phi=$PI*$ang/180;



    foreach(`cat CDs/$CDFile`){
	if(/^cell/){
	    print OUT;
	    chomp;
	    $CDheader=$_;
	    next;
	}
	chomp;
	my($CT, $cell, $time, $e1, $e2, $e3, $e4, $e5, $z, $x, $y, $d, $gweight) = split(/,/);
	unless($cell eq $prev){
	    my $standard = $division{$prev} - $birth{$prev};
	    if($standard=~/\d/ && $standard >0){
		my $observed = $prevTime - $birthTime;
#		print "$cell $observed $standard, " . $standard/$observed . "\n";
		if($observed>0){
		    push(@ratios, $standard/$observed);
		}
	    }
	    
	    $prev = $cell;
	    $birthTime = $time;
	    $thisMinTime{$cell}=$time;
	}
	$prevTime = $time;
	my $oldX = $x;

	#recenter to (0,0,0)
	$x=-$xc+$x;
	$y=-$yc+$y;
	$z=-$zc+$z;
	
	#derotate (CONSIDER DESTRETCIHNG AS WELL)

	my $newX= $x*cos($phi) - $y*sin($phi);
	my $newY= $y*cos($phi) + $x*sin($phi);
	$x=$newX;
	$y=$newY;

	

	#recenter to (356,256,33.5)
	$x=$x+356;
	$y=$y+256;
	$z=$z+33.5;

	#flip axis to match ADL standard
	if($axis eq "AVR"){
	    $z= 67-$z;
	    $y=512-$y;
	}elsif($axis eq "PVL"){
	    $x=712-$x;
	    $y=512-$y;
	}elsif($axis eq "PDR"){	    
	    $x=712-$x;
	    $z=67-$z;
	}elsif($axis eq "ADL"){
	}else{
	    print STDERR "Couldn't parse axis $axis\n";
	}

	#need to check the above for logic with different orientations
	# update - appears correct now.  May get better alignment with stretching - see below
	
#	print "$name $CT $axis $oldX $x $y $z \n";
#	print OUT join("\t", $CT, $cell, $time, $e1, $e2, $e3, $e4, $e5, $z, $x, $y, $d, $gweight) . "\n";

	

	#store in struture to allow time-warping
	if(exists($loc{$cell})){
	    $loc{$cell}= [@{$loc{$cell}}, [$e1, $e2, $e3, $e4, $e5, $z, $x, $y, $d, $gweight]];
	}else{
	    $loc{$cell} = [[$e1, $e2, $e3, $e4, $e5, $z, $x, $y, $d, $gweight]];
	}
	if($time > $thisMaxTime{$cell}){
	    $thisMaxTime{$cell}=$time;
	}

    }
    my $sum;
    foreach(@ratios){
	$sum += $_;
    }
    my $timeCorrection = $sum/int(@ratios);


    foreach(sort keys %loc){
	my $startOutIndex = $birth{$_};
	my $endOutIndex=$birth{$_};
	
	if($division{$_}=~/\d/){#Cell supposed to divide
	    my $length=$division{$_}-$birth{$_};
	    $endOutIndex=$division{$_};
	    if(!defined($loc{$daughter{$_}})){ #divider cell didn't actually divide in this series - fill partial or full trajectory based on rate
		my $altLength = int(1+$timeCorrection * ($thisMaxTime{$_}-$thisMinTime{$_}));
		if($altLength == 0){
		    $altLength=1;
		}
		if($altLength<$length){ 
		    $endOutIndex = $birth{$_} + $altLength;
		}
#		print STDERR "DIVISION NOT OBSERVED\n";
	    }else{
#		print STDERR "DIVIDED\n";
	    }
	}else{#cell not supposed to divide
	    my $altLength = int(1+$timeCorrection * ($thisMaxTime{$_}-$thisMinTime{$_}));
	    if ($altLength==0){
		$altLength=1;
	    }
	    $endOutIndex = $birth{$_} + $altLength;
#	    print STDERR "NOT SUPPOSED TO DIVIDE - UPDATED BASED ON LENGTH $altLength\n";
	}

	if(!defined($birth{$_})){#model doesn't know about this cell

	    #TODO: UPDATE MODEL TO KNOW ABOUT THESE CELLS!!!
#	    print STDERR "MODEL DOESN'T KNOW ABOUT $_ - SKIPPING\n";
	    next;
	}

	if(!defined($loc{$parent{$_}})){#birth not observed in this series - work backwards from end with rate
	    my $altLength = int(1+$timeCorrection * ($thisMaxTime{$_}-$thisMinTime{$_}));
	    if ($altLength==0){
		$altLength=1;
	    }
	    if($division{$_}-$altLength > $birth{$_}){
		$startOutIndex=$division{$_}-$altLength; # compress into shorter initial range
	    }
#	    print STDERR "BIRTH OF $_ NOT OBSERVED FOR $name - NEW RANGE $startOutIndex, $endOutIndex\n";
	}
	
#	    print STDERR "$_ didn't divide - interpolating\n";


	

	
	my $length=$endOutIndex-$startOutIndex;
#	print STDERR "$_ $length $startOutIndex $endOutIndex $birth{$_} $division{$_} THIS EMBRYO $name $thisMinTime{$_} $thisMaxTime{$_} \n";
	
	for(my $t=$startOutIndex; $t<=$endOutIndex; $t++){
	    my $maxIndex= $#{$loc{$_}};
	    my $thisIndex = int($maxIndex*($t-$startOutIndex)/$length);
#	    print STDERR "$t $maxIndex $thisIndex $length\n";
	    #OUTPUT TO ACD FILES FOR EACH SERIES TO ALLOW LATER CHECKING
	    print OUT join(",", "$_:$t", $_, $t,@{$loc{$_}->[$thisIndex]}) . "\n"; 
	    $data{"$_:$t"}{$name} = [$_, $t,@{$loc{$_}->[$thisIndex]}];
	}


    }
    close OUT;
}


#LOAD AND AVERAGE ACD (aligned CD) files
my $CT;
print $CDheader . "\n";
foreach $CT (sort ByCT keys %data){
#foreach $CT (sort keys %data){
#    print STDERR "$data{$CT}{$series[0]}[0] $data{$CT}{$series[0]}[1]\n"; 
    my $series;
    my ($cell, $time, $none,$global,$local,$blot,$cross,$z,$x,$y,$size,$gweight);

    my $count =0;
    foreach $series (sort keys %{$data{$CT}}){
	my @thisData = @{$data{$CT}{$series}};
	$cell =$thisData[0];
	$time = $thisData[1];
	$none+= $thisData[2];
	$global+= $thisData[3];
	$local+= $thisData[4];
	$blot+= $thisData[5];
	$cross+= $thisData[6];
	$z+= $thisData[7];
	$x+= $thisData[8];
	$y+= $thisData[9];
	$size+= $thisData[10];
	$gweight+= $thisData[11];

	$count++;
    }
#    print STDERR "$CT Z $z X $x Y $y D $size/$count GW $gweight/$count\n";

    print join(",", $CT, $cell, $time, int $none/$count, int $global/$count, int $local/$count, int $blot/$count, int $cross/$count, int(10*$z/$count)/10, int $x/$count, int $y/$count,int $size/$count,int $gweight/$count) . "\n";
}

sub ByCT{
#    my ($a, $b) = @_;
    my($aCell,$aTime) = split(/:/,$a);
    my($bCell,$bTime) = split(/:/,$b);
    if($aCell lt $bCell){
#	print STDERR "$a < $b\n";
	return -1;
    }elsif($aCell gt $bCell){
#	print STDERR "$a > $b\n";
	return 1;
    }else {
#	print STDERR "$a and $b same cell - sort by time\n";
	return $aTime<=> $bTime
    }
}
