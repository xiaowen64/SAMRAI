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

for t in ArrayDataNormOpsComplex ArrayDataNormOpsInteger \
HierarchyDataOpsComplex PatchSideDataOpsInteger \
PatchSideDataOpsComplex PatchNodeDataNormOpsComplex \
PatchSideDataNormOpsComplex HierarchySideDataOpsInteger \
HierarchySideDataOpsComplex HierarchyNodeDataOpsInteger \
HierarchyNodeDataOpsComplex HierarchyDataOpsManager \
HierarchyDataOpsInteger HierarchyDataOpsComplex \
PatchFaceDataOpsInteger PatchFaceDataOpsComplex \
PatchFaceDataNormOpsComplex HierarchyFaceDataOpsInteger \
PatchNodeDataOpsComplex HierarchyFaceDataOpsComplex \
PatchEdgeDataOpsInteger PatchEdgeDataOpsComplex \
PatchEdgeDataNormOpsComplex HierarchyEdgeDataOpsInteger \
HierarchyEdgeDataOpsComplex PatchCellDataOpsInteger \
PatchCellDataOpsComplex PatchCellDataNormOpsComplex \
HierarchyCellDataOpsInteger HierarchyCellDataOpsComplex \
PatchNodeDataOpsInteger 
do 
${MT} default.filenames ./tmp math $t NDIM 
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





