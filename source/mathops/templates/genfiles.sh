##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/templates/genfiles.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 2157 $
## Modified:    $LastChangedDate: 2008-04-28 13:55:06 -0700 (Mon, 28 Apr 2008) $
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

for t in ArrayDataNormOpsComplex \
HierarchyDataOpsComplex \
PatchSideDataOpsComplex \
PatchNodeDataNormOpsComplex \
PatchSideDataNormOpsComplex \
HierarchySideDataOpsComplex \
HierarchyNodeDataOpsComplex \
HierarchyDataOpsManager \
HierarchyDataOpsComplex \
PatchFaceDataOpsComplex \
PatchFaceDataNormOpsComplex \
PatchNodeDataOpsComplex \
HierarchyFaceDataOpsComplex \
PatchEdgeDataOpsComplex \
PatchEdgeDataNormOpsComplex \
HierarchyEdgeDataOpsComplex \
PatchCellDataOpsComplex \
PatchCellDataNormOpsComplex \
HierarchyCellDataOpsComplex 
do 
${MT} dcomplex.filenames ./tmp math $t NDIM 
done

for t in ArrayDataNormOpsInteger \
PatchSideDataOpsInteger \
HierarchySideDataOpsInteger \
HierarchyNodeDataOpsInteger \
HierarchyDataOpsManager \
HierarchyDataOpsInteger \
PatchFaceDataOpsInteger \
HierarchyFaceDataOpsInteger \
PatchEdgeDataOpsInteger \
HierarchyEdgeDataOpsInteger \
PatchCellDataOpsInteger \
HierarchyCellDataOpsInteger \
PatchNodeDataOpsInteger 
do 
${MT} int.filenames ./tmp math $t NDIM 
done

for t in int float double dcomplex; do
    ${MT} $t.filenames ./tmp math ArrayDataBasicOps NDIM\,$t
done

for t in double float; do
    ${MT} $t.filenames ./tmp math ArrayDataNormOpsReal NDIM\,$t
    ${MT} $t.filenames ./tmp math ArrayDataMiscellaneousOpsReal NDIM\,$t

    ${MT} $t.filenames ./tmp math HierarchyDataOpsReal NDIM\,$t
    ${MT} $t.filenames ./tmp tbox tbox::Pointer math::HierarchyDataOpsReal NDIM\,$t
    ${MT} $t.filenames ./tmp tbox tbox::Array tbox::Pointer math::HierarchyDataOpsReal NDIM\,$t
done

${MT} dcomplex.filenames ./tmp tbox tbox::Pointer math::HierarchyDataOpsComplex NDIM
${MT} default.filenames  ./tmp tbox tbox::Pointer math::HierarchyDataOpsInteger NDIM
${MT} dcomplex.filenames ./tmp tbox tbox::Array tbox::Pointer math::HierarchyDataOpsComplex NDIM
${MT} default.filenames  ./tmp tbox tbox::Array tbox::Pointer math::HierarchyDataOpsInteger NDIM
for g in Cell Edge Face Node Side; do 
    for t in int double float dcomplex; do
	${MT} $t.filenames ./tmp math math::Patch${g}DataBasicOps NDIM\,$t
    done

    for t in double float; do
	${MT} $t.filenames ./tmp math math::Patch${g}DataNormOpsReal NDIM\,$t
	${MT} $t.filenames ./tmp math math::Patch${g}DataMiscellaneousOpsReal NDIM\,$t
	${MT} $t.filenames ./tmp math math::Patch${g}DataOpsReal NDIM\,$t
	${MT} $t.filenames ./tmp tbox tbox::Pointer math::Patch${g}DataOpsReal NDIM\,$t
    done
    ${MT} dcomplex.filenames ./tmp tbox tbox::Pointer math::Patch${g}DataOpsComplex NDIM
    ${MT} default.filenames  ./tmp tbox tbox::Pointer math::Patch${g}DataOpsInteger NDIM
done
for g in Cell Edge Face Node Side; do
    for t in double float; do
	${MT} $t.filenames ./tmp math math::Hierarchy${g}DataOpsReal NDIM\,$t
	${MT} $t.filenames ./tmp tbox tbox::Pointer math::Hierarchy${g}DataOpsReal NDIM\,$t
    done 
    ${MT} dcomplex.filenames ./tmp tbox tbox::Pointer math::Hierarchy${g}DataOpsComplex NDIM
    ${MT} int.filenames      ./tmp tbox tbox::Pointer math::Hierarchy${g}DataOpsInteger NDIM
done


#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticXd ./tmp/*.C

sh ../../scripts/object.sh ./tmp automaticXd
sh ../../scripts/depend

rm -rf ./tmp


exit 0





