##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/mesh/templates/genfiles.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 1704 $
## Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
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
AsyncBergerRigoutsos AsyncBergerRigoutsosNode BalanceUtilities \
BaseGriddingAlgorithm BergerRigoutsos BoxGeneratorStrategy GriddingAlgorithm \
HistogramBox LoadBalancer LoadBalanceStrategy MultiblockGriddingAlgorithm \
MultiblockGriddingTagger StandardTagAndInitialize StandardTagAndInitStrategy \
TagAndInitializeStrategy 
do
${MT} default.filenames ./tmpXd mesh $t NDIM
done

# NonDim
${MT} default.filenames ./tmpNond tbox tbox::Array mesh::SpatialKey

${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::BaseGriddingAlgorithm NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::BergerRigoutsos NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::BoxGeneratorStrategy NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::GriddingAlgorithm NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::LoadBalancer NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::LoadBalanceStrategy NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::StandardTagAndInitialize NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::StandardTagAndInitStrategy NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::TagAndInitializeStrategy NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::AsyncBergerRigoutsos NDIM
${MT} default.filenames ./tmpXd tbox tbox::Pointer mesh::MultiblockGriddingAlgorithm NDIM

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





