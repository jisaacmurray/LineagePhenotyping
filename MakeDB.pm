#!/usr/bin/perl 
use strict;



#DB goes to hash describing the xml file structure.  Standard keys include:

sub MakeDB{
    my $dbDir = "/gpfs/fs0/l/murr/embryoDB";
    my %DB;

    my %constructs;

    print("Reading strains\n");
    foreach(`cat /gpfs/fs0/l/murr/Strains.csv`){
	chomp;
	$_ =~tr/\"//;
	my ($strain, $keep, $reporter, $allele,$construct, @etc) = split(/,/);
	$constructs{$strain} = $construct;
#	print STDERR;
    }



    print("Reading db entries\n");

    foreach(`ls $dbDir`){
#	print;
	chomp;
	
	if(/xml$/){
	    my %thisEntry;
	    open(FILE, $dbDir . "/" . $_);
	    foreach(<FILE>){
		if(/<(.+)\s.+=\"(.+)\"\/>/){
		    $thisEntry{$1} = $2;
		}
	    }
	    if($thisEntry{"comments"} =~/ENCODE/){
		$thisEntry{"type"} = "protein";
	    }else{
		$thisEntry{"type"} = "promoter";
	    }
	    if(exists($constructs{$thisEntry{"strain"}})){
		$thisEntry{"construct"} = $constructs{$thisEntry{"strain"}};
	    }else{
		$thisEntry{"construct"} = "NULL";
	    }


	    if($thisEntry{"redsig"} =~/hnd/){
#		print join (" " , $_ , $thisEntry{"series"}, $thisEntry{"redsig"}) . "\n";
	    }

	    my $auxInfo = $thisEntry{"annots"} . "/dats/" . $thisEntry{"series"} . "AuxInfo.csv";
	    my $line=0;
	    my @keys;
	    my @values;
	    if (-e $auxInfo){
		foreach(`cat $auxInfo`){
		    chomp;
		    if ($line==0){
			@keys = split(/,/);
			$line++;
		    }else{
			@values = split(/,/);
		    }
		}
		
		for(my $i=0;$i<=$#keys;$i++){
		    if($keys[$i] eq "axis"){
			$thisEntry{$keys[$i]} = uc($values[$i]);
			
		    }else{
			$thisEntry{$keys[$i]} = $values[$i];
		    }
#		    print STDERR join(" ", $thisEntry{"series"}, $keys[$i], $values[$i]) . "\n";
		}
		
		
	    }
	    
#	    my $CA =  $thisEntry{"annots"} . "/dats/CA" . $thisEntry{"series"} . ".csv";
	    $thisEntry{"ABala"} = 0;
	    $thisEntry{"C"} = 0;
# 	    if(-e $CA){
# 		foreach(`cat $CA`){
# 		    my @line = split(/,/);
# 		    if($line[1] eq "C"){
# #			print STDERR join(" ", $thisEntry{"series"}, "ABala", $line[2]) . "\n";
			
# 			$thisEntry{"C"} = $line[2];
# 			last;
# 		    }else{
# #			print STDERR "$line[1] not ABala for" . $thisEntry{"series"}. "\n";
# 		    }
# 		}
# 	    }
# #	    print join("\t", $thisEntry{"series"}, $thisEntry{"axis"}) . "\n";
    
	    $DB{$thisEntry{"series"}} = \%thisEntry;
	}
    }
    return \%DB
}


sub IsEdited{
    my ($DBref, $series) = @_;
    if(
	($DBref->{$series}{"status"} ne "deleted")&&
#	(($DBref->{$series}{"editedby"} =~/JIM/)|| ($DBref->{$series}{"checkedby"} =~ "JIM")) &&
	($DBref->{$series}{"editedtimepts"} =~ /^\d+$/) && #require edited time be a number
#	($DBref->{$series}{"redsig"} ne "none") && # filter out no-reporter series
	($DBref->{$series}{"editedcells"} >= 40) &&  #require edited cells be 40 or more
	($DBref->{$series}{"date"} >= 20100101) && #require generated at murray lab
	($DBref->{$series}{"redsig"} =~/^[\w\-\.]+$/) ##legal gene name 
       
	){
	if(($DBref->{$series}{"construct"} =~/\?/) || ($DBref->{$series}{"redsig"} =~/\?/) || ($DBref->{$series}{"construct"} =~/pJIM11/)){
#	    print STDERR "SKIPPING SERIES $series WITH ambiguous construct " . $DBref->{$series}{"construct"} . "\n";
	    return 0;
	}
	else{
	    return 1;
	}
    } else {
#	print STDERR "$series Not Edited!\n";
	return 0;
    }
}



sub IsStandard{
    
    my ($DBref, $series) = @_;

    if(
	($DBref->{$series}{"status"} ne "deleted")&&
#	(($DBref->{$series}{"editedby"} =~/JIM/)|| ($DBref->{$series}{"checkedby"} =~ "JIM")) &&
	($DBref->{$series}{"redsig"} ne "none") && # filter out no-reporter series
	($DBref->{$series}{"editedtimepts"} =~ /^\d+$/) && #require edited time be a number
	($DBref->{$series}{"editedcells"} >= 40) && #require edited cells be 40 or more
#	($DBref->{$series}{"treatments"} eq "none" || $DBref->{$series}{"treatments"} eq 'n/a')  &&#remove RNAi series
	($DBref->{$series}{"date"} >= 20100101) && #require generated at murray lab
	($DBref->{$series}{"redsig"} =~/^[\w\-\.]+$/) ##legal gene name 
	){
	if(($DBref->{$series}{"construct"} =~/\?/) || ($DBref->{$series}{"redsig"} =~/\?/) || ($DBref->{$series}{"construct"} =~/pJIM11/)){
#	    print STDERR "SKIPPING SERIES $series WITH ambiguous or bad construct " . $DBref->{$series}{"construct"} . "\n";
	    return 0;
	}
	else{
	    return 1;
	}
	
    } else {
	return 0;
    }
}

return 1;

#test code
# my %DB = %{MakeDB()};

# foreach(sort {$DB{$a}{"redsig"} cmp $DB{$b}{"redsig"}} keys (%DB)){

#     if($DB{$_}{"editedcells"} =~/^\d+$/ && 
#        $DB{$_}{"editedcells"}>100 &&
#        $DB{$_}{"redsig"} ne "none" &&
#        $DB{$_}{"person"} eq "murray" &&
#        $DB{$_}{"redsig"} ne "n\\a" ){
# 	print join("\t", 
# 		   $DB{$_}{"redsig"}, 
# 		   $DB{$_}{"strain"}, 
# 		   $DB{$_}{"treatments"}, 
# 		   $DB{$_}{"person"}, 
# 		   $DB{$_}{"editedtimepts"},
# 		   $DB{$_}{"editedcells"},
# 		   $_, 
# 		   $DB{$_}{"comments"}
		   

# 	    ) . "\n";
#     }
# }
