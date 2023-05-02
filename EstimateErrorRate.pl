#!/usr/bin/perl -w
use MakeDB;
use strict;

my $DBref = MakeDB();
my @series;
foreach(<>){
    chomp;
    push(@series, $_);
}

print join ("\t",qw(series file tp fp fn move)) . "\n";
my $series;
foreach $series (@series){
    my $annots = $DBref->{$series}{"annots"};
    my $zip = $annots . "/dats/" . $series . ".zip";
    my $editzip = $annots . "/dats/" . $series . "-edit.zip";

    if(!-e $zip || !-e $editzip){
	print STDERR"Can't get files for $series\n";
	next;
    }

    my $editedTP = $DBref->{$series}{"editedtimepts"};
    my $checked = $DBref->{$series}{"checkedby"};
    if($checked=~/(\d+),/){
	$editedTP=$1;
    }
#    print join ("\t", $series, $editedTP, $zip) . "\n";

    system "rm -Rf TMPzip";
    system "rm -Rf TMPeditzip";
    mkdir("TMPzip");
    mkdir("TMPeditzip");
    chdir("TMPzip");
#    print `pwd`;
    if(system("unzip -qq $zip 2>TMPOUT")){
#	print "CAN'T OPEN ZIPFILE $zip FOR $series - SKIPPING\n";
	chdir "..";
	next;
    }
    chdir("../TMPeditzip");
#    print `pwd`;
    if(system("unzip -qq $editzip 2>TMPOUT")){
#	print  "CAN'T OPEN EDITZIPFILE $editzip FOR $series - SKIPPING\n";
	chdir "..";
	next;
    }
    chdir("..");
    print  STDERR "SUCCESSFULLY EXTRACTED FILES FOR $series\n";

    my $nucFile;
    foreach $nucFile (`ls TMPzip/nuclei`){
	chomp $nucFile;
	if($nucFile=~/t(\d\d\d)-nuclei/){
	    if($1 <= $editedTP){
		my $c1=0;
		my $c2=0;
		my %pos;
		my %posEdit;
		my %names;
		my %names2;
		foreach(`cat TMPzip/nuclei/$nucFile`){
		    $c1++;
		    my ($index,$status,$parent, $d1,$d2,$x,$y,$z,$d,$name,@etc) = split(/,\s+/);
		    $index = int($index);
		    if($status==1){
			$pos{$index}=[$x,$y,$z];
			$names{$index}=$name;
#			print "ORIG $index,x\n";
		    }
		}
		foreach(`cat TMPeditzip/nuclei/$nucFile`){
		    $c2++;
		    my ($index,$status,$parent, $d1,$d2,$x,$y,$z,$d,$name,@etc) = split(/,\s+/);
		    $index = int($index);
		    if($status==1){
			$posEdit{$index} = [$x,$y,$z];
#			print "EDIT $index,x\n";
			$names2{$index}=$name;
		    }
		}
		
		my $fp=0;
		my $tp=0;
		my $fn=0;
		my $move=0;

		foreach(sort keys %pos){
		    if(defined($posEdit{$_})){
			$tp++;
			if($posEdit{$_}[0] != $pos{$_}[0] ||
			   $posEdit{$_}[1] != $pos{$_}[1] ||
			   $posEdit{$_}[2] != $pos{$_}[2]){
			    $move++;
			}
		    }else{
			$fp++;
		    }
		}
		foreach(sort keys %posEdit){
		    if(!defined($pos{$_})){
			$fn++;
		    }
		}
		print join ("\t",$series, $nucFile, $tp, $fp, $fn,$move) . "\n";
	    }
	}
    }
#    compare nuc file by nuc file
}
