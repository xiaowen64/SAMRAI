#!/usr/bin/perl
##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/inputdb/scanner_fixup.pl $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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

    s/#if YY_STACK_USED/#ifdef YY_STACK_USED/;
    s/#if YY_ALWAYS_INTERACTIVE/#ifdef YY_ALWAYS_INTERACTIVE/;
    s/#if YY_MAIN/#ifdef YY_MAIN/;

    print $_;
}
