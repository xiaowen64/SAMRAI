#!/usr/local/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/scripts/cmp.pl $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
## Description: perl script to compare two files but ignore CVS comments
##

##
## Usage: cmp.pl <file1> <file2>
##

$ANAME = shift(@ARGV);
$BNAME = shift(@ARGV);

open(AFILE, "$ANAME") || die "Cannot open input file $ANAME...";
open(BFILE, "$BNAME") || die "Cannot open input file $BNAME...";

while (!eof(AFILE) && !eof(BFILE)) {
   $ALINE = <AFILE>;
   $BLINE = <BFILE>;
   $_ = $ALINE;

   if (!/^(\/\/|c|C|#|##| \*)[ ]*(Release:[\t ]*\$Name|Revision:[\t ]*\$LastChangedRevision|Modified:[\t ]*\$LastChangedDate):[^\$]*\$/o) {
      if ($ALINE ne $BLINE) {
         exit 1;
      }
   }
}

exit 0 if (eof(AFILE) && eof(BFILE));
exit 1;
