##
## File:        genfiles.sh
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2005 The Regents of the University of California
## Revision:    $Revision$
## Modified:    $Date$
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

for t in CoarsenAlgorithm CoarsenClasses CoarsenCopyTransaction \
    CoarsenPatchStrategy CoarsenSchedule CoarsenOperator \
    CoarsenTransactionFactory StandardCoarsenTransactionFactory   \
    RefineAlgorithm RefineClasses RefineCopyTransaction \
    RefinePatchStrategy  RefineSchedule RefineTimeTransaction   \
    RefineTransactionFactory StandardRefineTransactionFactory \
    RefineOperator TimeInterpolateOperator FillBoxSet \
    Geometry \
    FillBoxSet LocallyActiveDataFillBox LocallyActiveDataFillBoxSet \
    LocallyActiveDataCoarsenPatchStrategy LocallyActiveDataCoarsenAlgorithm \
    LocallyActiveDataCoarsenTransactionFactory \
    StandardLocallyActiveDataCoarsenTransactionFactory LocallyActiveDataCoarsenSchedule \
    LocallyActiveDataRefinePatchStrategy LocallyActiveDataRefineAlgorithm \
    LocallyActiveDataRefineTransactionFactory \
    StandardLocallyActiveDataRefineTransactionFactory LocallyActiveDataRefineSchedule
do
  ${MT} default.filenames ./tmp xfer $t NDIM
done


${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Array tbox::Pointer xfer::RefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Pointer xfer::RefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array xfer::FillBoxSet NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Pointer xfer::LocallyActiveDataRefineSchedule NDIM

${MT} default.filenames ./tmp tbox tbox::Array xfer::FillBoxSet NDIM

${MT} default.filenames ./tmp tbox tbox::Array xfer::LocallyActiveDataFillBoxSet NDIM

${MT} default.filenames ./tmp tbox tbox::Array tbox::List xfer::CoarsenClasses\<NDIM\>::Data
${MT} default.filenames ./tmp tbox tbox::Array tbox::List xfer::RefineClasses\<NDIM\>::Data

${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer xfer::CoarsenSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer xfer::LocallyActiveDataCoarsenSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer xfer::RefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer xfer::LocallyActiveDataRefineSchedule NDIM

${MT} default.filenames ./tmp tbox tbox::List xfer::LocallyActiveDataFillBox NDIM 

${MT} default.filenames ./tmp tbox tbox::List xfer::CoarsenClasses\<NDIM\>::Data
${MT} default.filenames ./tmp tbox tbox::List xfer::RefineClasses\<NDIM\>::Data

${MT} default.filenames ./tmp tbox tbox::List tbox::Pointer xfer::CoarsenOperator NDIM
${MT} default.filenames ./tmp tbox tbox::List tbox::Pointer xfer::RefineOperator NDIM
${MT} default.filenames ./tmp tbox tbox::List tbox::Pointer xfer::RefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::List tbox::Pointer xfer::LocallyActiveDataRefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::List tbox::Pointer xfer::TimeInterpolateOperator NDIM


${MT} default.filenames ./tmp tbox tbox::Pointer xfer::Geometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::CoarsenAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::CoarsenClasses NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::CoarsenOperator NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::CoarsenSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::CoarsenTransactionFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::StandardCoarsenTransactionFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::RefineAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::RefineClasses NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::RefineOperator NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::RefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::RefineTransactionFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::StandardRefineTransactionFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::TimeInterpolateOperator NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::LocallyActiveDataCoarsenAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::LocallyActiveDataCoarsenSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::LocallyActiveDataCoarsenTransactionFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::StandardLocallyActiveDataCoarsenTransactionFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::LocallyActiveDataRefineAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::LocallyActiveDataRefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::LocallyActiveDataRefineTransactionFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::StandardLocallyActiveDataRefineTransactionFactory NDIM

#
# other list templates
#

#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticXd ./tmp/*.C

sh ../../scripts/object.sh ./tmp automaticXd
sh ../../scripts/depend

rm -rf ./tmp


exit 0

