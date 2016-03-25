#!/usr/bin/perl

while(<>) {
    s/^#line.*//;
    s/.*Revision:.*//;
    s/.*Date:.*//;
    s/.*Header:.*//;

    s/YY_DO_BEFORE_ACTION;/YY_DO_BEFORE_ACTION/;
    s/^(\s)+;$/$1do {} while(0);/;

    print $_;
}
