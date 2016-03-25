##
## File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/inputdb/genfiles.sh $
## Package:	SAMRAI toolbox
## Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:	$LastChangedRevision: 1917 $
## Modified:	$LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
## Description:	simple shell script to generate flex and bison files
##

dir_name=`echo ${0} | sed -e 's:^\([^/]*\)$:./\1:' -e 's:/[^/]*$::'`;
cd $dir_name

#
# Use yacc since ASCI red does not support alloca() function used by bison
#

bison -d -p yy Grammar.y
perl grammer_fixup.pl Grammar.tab.c > Grammar.C
perl grammer_fixup.pl Grammar.tab.h > Grammar.h
rm Grammar.tab.c
rm Grammar.tab.h

#
# Scanner requires flex due to input reading method with MPI
#

flex -Pyy -otemp.$$ Scanner.l
perl scanner_fixup.pl temp.$$ > Scanner.C 
rm temp.$$
