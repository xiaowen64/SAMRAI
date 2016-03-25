#!/usr/apps/bin/perl

## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/scripts/source_manipulation/removeNameKeyword.pl $
## Package:     SAMRAI application
## Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Script in source-manipulation utility.


# This script removes the <DIM> from constructor method names

require "find.pl";
use Cwd;
use File::Basename;
use Fcntl ':mode';

&find('.');

exit;

sub wanted {
    if ($File::Find::dir =~ m/.svn/o ) {
	$File::Find::prune = true;
    } elsif ( m/.svn/ ) {
    } elsif ( -f ) {
	doit($name);
    }
}

sub doit {
    my ($filename) = shift;
    print("$filename\n") if $debug;

    # The doit is invoked in the subdir so need to strip off 
    # the directory
    $filename=basename $filename;

    $ofilename=$filename . ".tmp";

    $mode = S_IMODE((stat($filename))[2]);

    stat($filename);
    if( -x _ ) {
	$executable=1;
    } else{
	$executable=0;
    }

    open(IN, $filename) || die "Cannot open input file $filename : $!";
    open OUT, "> $ofilename" || die "Cannot open output file $ofilename";

     while (<IN>) {
	 } else {
	     print OUT;
	 }
     }

     close IN  || die "Cannot close file $filename";
     close OUT || die "Cannot close file $ofilename";

    unlink($filename);
    rename($ofilename, $filename);

    chmod $mode, $filename;

}
