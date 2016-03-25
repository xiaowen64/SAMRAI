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

#
# The templates for the NDIM classes in this package
#

for t in \
MBUtilities MultiblockPatchLevel MultiblockPatchHierarchy \
MultiblockGriddingAlgorithm MultiblockGriddingTagger \
MultiblockCoarsenAlgorithm MultiblockCoarsenSchedule MultiblockCoarsenPatchStrategy \
MultiblockRefineAlgorithm MultiblockRefineSchedule MultiblockRefinePatchStrategy

do 
${MT} default.filenames ./tmp mblk $t NDIM 
done

${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Pointer mblk::MultiblockPatchLevel NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Pointer mblk::MultiblockRefineSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::List mblk::MultiblockPatchHierarchy\<NDIM\>::Neighbor
${MT} default.filenames ./tmp tbox tbox::Array tbox::List mblk::MultiblockRefineSchedule\<NDIM\>::SingularityPatch
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer mblk::MultiblockPatchLevel NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer mblk::MultiblockRefineSchedule NDIM

${MT} default.filenames ./tmp tbox tbox::List mblk::MultiblockPatchHierarchy\<NDIM\>::Neighbor
${MT} default.filenames ./tmp tbox tbox::List mblk::MultiblockRefineSchedule\<NDIM\>::SingularityPatch

${MT} default.filenames ./tmp tbox tbox::Pointer mblk::MultiblockCoarsenAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer mblk::MultiblockCoarsenSchedule NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer mblk::MultiblockGriddingAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer mblk::MultiblockPatchHierarchy NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer mblk::MultiblockPatchLevel NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer mblk::MultiblockRefineAlgorithm NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer mblk::MultiblockRefineSchedule NDIM

for t in double float bool dcomplex int; do
    ${MT} $t.filenames ./tmp mblk mblk::MBDataUtilities NDIM\,$t
done


#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticXd ./tmp/*.C

sh ../../scripts/object.sh ./tmp automaticXd
sh ../../scripts/depend

rm -rf ./tmp


exit 0





