##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/transfer/templates/genfiles.sh $
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

for t in CoarsenAlgorithm CoarsenClasses CoarsenCopyTransaction \
    CoarsenOperator CoarsenPatchStrategy CoarsenSchedule \
    CoarsenTransactionFactory FillBoxSet Geometry \
    LocallyActiveDataCoarsenPatchStrategy LocallyActiveDataCoarsenAlgorithm \
    LocallyActiveDataCoarsenSchedule \
    LocallyActiveDataCoarsenTransactionFactory LocallyActiveDataFillBox \
    LocallyActiveDataFillBoxSet LocallyActiveDataRefinePatchStrategy \
    LocallyActiveDataRefineAlgorithm LocallyActiveDataRefineSchedule \
    LocallyActiveDataRefineTransactionFactory MultiblockCoarsenAlgorithm \
    MultiblockCoarsenPatchStrategy MultiblockCoarsenSchedule \
    MultiblockRefineAlgorithm MultiblockRefinePatchStrategy \
    MultiblockRefineSchedule RefineAlgorithm RefineClasses \
    RefineCopyTransaction RefineOperator RefinePatchStrategy RefineSchedule \
    RefineTimeTransaction RefineTransactionFactory \
    StandardCoarsenTransactionFactory \
    StandardLocallyActiveDataCoarsenTransactionFactory \
    StandardLocallyActiveDataRefineTransactionFactory \
    StandardRefineTransactionFactory TimeInterpolateOperator
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

${MT} default.filenames ./tmp tbox tbox::Pointer xfer::MultiblockCoarsenAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::MultiblockCoarsenSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Pointer xfer::MultiblockRefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::List xfer::MultiblockRefineSchedule\<NDIM\>::SingularityPatch
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer xfer::MultiblockRefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::List xfer::MultiblockRefineSchedule\<NDIM\>::SingularityPatch
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::MultiblockRefineAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer xfer::MultiblockRefineSchedule NDIM


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

