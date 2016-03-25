##
## File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/inputdb/genfiles.sh $
## Package:	SAMRAI toolbox
## Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:	$LastChangedRevision: 1704 $
## Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
