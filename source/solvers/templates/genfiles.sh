##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/templates/genfiles.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1768 $
## Modified:    $LastChangedDate: 2007-12-11 16:02:04 -0800 (Tue, 11 Dec 2007) $
## Description: shell script to create SAMRAI template files in the repository
##

dir_name=`echo ${0} | sed -e 's:^\([^/]*\)$:./\1:' -e 's:/[^/]*$::'`;
cd $dir_name

#
# set up PERL program name and the make-templates.pl command
#

PERL=${PERL:-perl}
MT="$PERL ../../scripts/make-template.pl"

#
# create a new scratch temporary directory
#

rm -rf ./tmpXd
mkdir ./tmpXd

rm -rf ./tmpNond
mkdir ./tmpNond

#
# Create empty filename files for each type
#
for t in default Complex Double Float Integer bool dcomplex char float double int all
do
    touch ./tmpXd/$t.filenames
    touch ./tmpNond/$t.filenames
done

#
# The templates for the NDIM classes in this package
#

for t in \
FACOperatorStrategy  FACPreconditioner \
NonlinearSolverStrategy SNES_SAMRAIContext KINSOL_SAMRAIContext \
CartesianRobinBcHelper CellPoissonFACOps CellPoissonFACSolver \
CellPoissonHypreSolver GhostCellRobinBcCoefs LocationIndexRobinBcCoefs \
RobinBcCoefStrategy SimpleCellRobinBcCoefs \
Sundials_SAMRAIVector
do
${MT} default.filenames ./tmpXd solv $t NDIM
done

for t in double float; do
    ${MT} $t.filenames ./tmpXd solv SAMRAIVectorReal NDIM\,$t
    ${MT} $t.filenames ./tmpXd tbox tbox::Pointer solv::SAMRAIVectorReal NDIM\,$t
done
for t in double; do
    ${MT} $t.filenames ./tmpXd solv PETSc_SAMRAIVectorReal NDIM\,$t 
done

for t in double; do
   ${MT} $t.filenames ./tmpNond solv PETScAbstractVectorReal $t
done

#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticXd ./tmpXd/*.C
sh ../../scripts/copy-if-change ./automaticNond ./tmpNond/*.C

sh ../../scripts/object.sh ./tmpXd automaticXd
sh ../../scripts/object.sh ./tmpNond automaticNond
sh ../../scripts/depend

rm -rf ./tmpXd
rm -rf ./tmpNond

exit 0





