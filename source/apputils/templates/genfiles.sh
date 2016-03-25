##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/apputils/templates/genfiles.sh $
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
CartesianVizamraiDataWriter  VisItDataWriter \
VisDerivedDataStrategy       VisMaterialsDataStrategy \
CubesPatchInterface          ElevenPatchInterface \
BoundaryNode                 CutCell \
EmbeddedBoundaryGeometry     \
EmbeddedBoundaryShape        \
EmbeddedBoundaryShapePolygon \
EmbeddedBoundaryShapeSphere
do 
${MT} default.filenames ./tmp appu $t NDIM 
done

${MT} default.filenames ./tmp tbox tbox::List appu::CartesianVizamraiDataWriter\<NDIM\>::VizamraiItem NDIM
${MT} default.filenames ./tmp tbox tbox::List appu::VisItDataWriter\<NDIM\>::VisItItem
${MT} default.filenames ./tmp tbox tbox::Pointer appu::CartesianVizamraiDataWriter NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer appu::VisItDataWriter NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer appu::EmbeddedBoundaryGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer appu::EmbeddedBoundaryShape NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer appu::CubesPatchInterface NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer appu::ElevenPatchInterface NDIM
${MT} default.filenames ./tmp tbox tbox::Array appu::BoundaryNode NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer appu::EmbeddedBoundaryShape NDIM

#
# Templates for CutCell
#

${MT} default.filenames ./tmp pdat pdat::IndexData NDIM\,appu::CutCell\<NDIM\>
${MT} default.filenames ./tmp pdat pdat::IndexDataFactory NDIM\,appu::CutCell\<NDIM\>
${MT} default.filenames ./tmp pdat pdat::IndexVariable NDIM\,appu::CutCell\<NDIM\>
${MT} default.filenames ./tmp tbox tbox::Array appu::CutCell\<NDIM\>
${MT} default.filenames ./tmp tbox tbox::Pointer appu::CutCell\<NDIM\>
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer appu::CutCell\<NDIM\>
${MT} default.filenames ./tmp tbox tbox::Pointer pdat::IndexData NDIM\,appu::CutCell\<NDIM\>
${MT} default.filenames ./tmp tbox tbox::Pointer pdat::IndexVariable NDIM\,appu::CutCell\<NDIM\>
${MT} default.filenames ./tmp tbox tbox::Pointer pdat::IndexDataFactory NDIM\,appu::CutCell\<NDIM\>


#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticXd ./tmp/*.C

sh ../../scripts/object.sh ./tmp automaticXd
sh ../../scripts/depend

rm -rf ./tmp


exit 0





