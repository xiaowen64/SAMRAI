##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/algorithm/templates/genfiles.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1917 $
## Modified:    $LastChangedDate: 2008-01-25 13:28:01 -0800 (Fri, 25 Jan 2008) $
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

rm -rf ./tmp
mkdir ./tmp

#
# Create empty filename files for each type
#
for t in default Complex Double Float Integer bool dcomplex char float double int all
do
    touch ./tmp/$t.filenames
done

#
# The templates for the NDIM classes in this package
#

for t in \
TimeRefinementIntegrator TimeRefinementLevelStrategy \
HyperbolicPatchStrategy HyperbolicLevelIntegrator \
ImplicitEquationStrategy  ImplicitIntegrator \
MblkPatchBoundaryNodeSum \
MethodOfLinesIntegrator MethodOfLinesPatchStrategy \
OuternodeSumTransaction OuternodeSumTransactionFactory  PatchBoundaryNodeSum \
OuteredgeSumTransaction OuteredgeSumTransactionFactory  PatchBoundaryEdgeSum \
LocallyActiveDataOuteredgeSumTransactionFactory \
LocallyActiveDataOuternodeSumTransactionFactory \
LocallyActiveDataPatchBoundaryEdgeSum LocallyActiveDataPatchBoundaryNodeSum
do
${MT} default.filenames ./tmp algs $t NDIM
done

${MT} default.filenames ./tmp tbox tbox::Pointer algs::HyperbolicLevelIntegrator NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer algs::MethodOfLinesIntegrator NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer algs::TimeRefinementIntegrator NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer algs::TimeRefinementLevelStrategy NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer algs::MblkPatchBoundaryNodeSum NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer algs::PatchBoundaryNodeSum NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer algs::PatchBoundaryEdgeSum NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer algs::PatchBoundaryEdgeSum NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer algs::PatchBoundaryNodeSum NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer algs::LocallyActiveDataPatchBoundaryNodeSum NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer algs::LocallyActiveDataPatchBoundaryEdgeSum NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer algs::LocallyActiveDataPatchBoundaryEdgeSum NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer algs::LocallyActiveDataPatchBoundaryNodeSum NDIM

#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticXd ./tmp/*.C

sh ../../scripts/object.sh ./tmp automaticXd
sh ../../scripts/depend

rm -rf ./tmp


exit 0





