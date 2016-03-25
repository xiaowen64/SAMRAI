#!/usr/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/inputdb/scanner_fixup.pl $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description: Script in input database package.
##

while(<>) {
    s/^#line.*//;
    s/.*Revision:.*//;
    s/.*Date:.*//;
    s/.*Header:.*//;

    # substitution to replace [yylval] with SAMRAI_[yylval]
    s/yylval/SAMRAI_yylval/g;

    s/YY_DO_BEFORE_ACTION;/YY_DO_BEFORE_ACTION/;
    s/^(\s)+;$/$1do {} while(0);/;

    print $_;
}
