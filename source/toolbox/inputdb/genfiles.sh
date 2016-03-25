##
## File:	make-parser.sh
## Package:	SAMRAI toolbox
## Copyright:	(c) 1997-2005 The Regents of the University of California
## Revision:	$Revision$
## Modified:	$Date$
## Description:	simple shell script to generate flex and bison files
##

dir_name=`echo ${0} | sed -e 's:^\([^/]*\)$:./\1:' -e 's:/[^/]*$::'`;
cd $dir_name

#
# Use yacc since ASCI red does not support alloca() function used by bison
#

bison -d -p yy Grammar.y
perl grammer_fixup.pl Grammar.tab.c > Grammar.C
rm Grammar.tab.c
mv Grammar.tab.h  Grammar.h

#
# Scanner requires flex due to input reading method with MPI
#

flex -Pyy -otemp.$$ Scanner.l
perl scanner_fixup.pl temp.$$ > Scanner.C 
rm temp.$$
