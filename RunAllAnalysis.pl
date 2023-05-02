#!/usr/bin/perl -w
use strict;

foreach(@ARGV){
    system "./CompareDivTime.pl $_";
    system "./ComparePositions.pl $_";
#Not fully implemented yet
#    system "./CompareMigration.pl $_";
    system "./GetAngles_revRotate.pl $_";
}
