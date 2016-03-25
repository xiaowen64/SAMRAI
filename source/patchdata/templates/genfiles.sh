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
ArrayDataIterator \
CellIndex CellIterator \
CellGeometry CellOverlap \
EdgeGeometry EdgeOverlap \
FaceGeometry FaceOverlap \
NodeGeometry NodeOverlap \
OuteredgeGeometry \
OuterfaceGeometry \
OuternodeGeometry \
OutersideGeometry \
SideGeometry SideOverlap \
EdgeIndex EdgeIterator \
SideIndex SideIterator \
NodeIndex NodeIterator  \
FaceIndex FaceIterator \
CellComplexLinearTimeInterpolateOp \
CellDoubleLinearTimeInterpolateOp CellFloatLinearTimeInterpolateOp \
OutersideComplexLinearTimeInterpolateOp \
OutersideDoubleLinearTimeInterpolateOp \
OutersideFloatLinearTimeInterpolateOp \
EdgeComplexLinearTimeInterpolateOp EdgeDoubleLinearTimeInterpolateOp \
EdgeFloatLinearTimeInterpolateOp SideComplexLinearTimeInterpolateOp \
SideDoubleLinearTimeInterpolateOp SideFloatLinearTimeInterpolateOp \
NodeComplexLinearTimeInterpolateOp NodeDoubleLinearTimeInterpolateOp \
NodeFloatLinearTimeInterpolateOp \
OuterfaceComplexLinearTimeInterpolateOp \
OuterfaceDoubleLinearTimeInterpolateOp \
OuterfaceFloatLinearTimeInterpolateOp \
FaceComplexLinearTimeInterpolateOp FaceDoubleLinearTimeInterpolateOp \
FaceFloatLinearTimeInterpolateOp CellComplexConstantRefine \
CellDoubleConstantRefine CellFloatConstantRefine \
CellIntegerConstantRefine SideComplexConstantRefine \
SideDoubleConstantRefine SideFloatConstantRefine \
SideIntegerConstantRefine EdgeComplexConstantRefine \
EdgeDoubleConstantRefine EdgeFloatConstantRefine \
EdgeIntegerConstantRefine OuternodeDoubleConstantCoarsen \
NodeComplexConstantAverage \
NodeDoubleConstantAverage NodeFloatConstantAverage \
NodeIntegerConstantAverage OuterfaceComplexConstantRefine \
OuterfaceDoubleConstantRefine OuterfaceFloatConstantRefine \
OuterfaceIntegerConstantRefine FaceComplexConstantRefine \
FaceDoubleConstantRefine FaceFloatConstantRefine \
FaceIntegerConstantRefine
do
  ${MT} default.filenames ./tmp pdat $t NDIM
done

for g in Node; do
    for t in Complex Double Float Integer; do
	${MT} $t.filenames ./tmp tbox tbox::Pointer pdat::${g}${t}ConstantAverage NDIM
    done
done
for g in Outernode; do
    for t in Double; do
	${MT} $t.filenames ./tmp tbox tbox::Pointer pdat::${g}${t}ConstantCoarsen NDIM
    done
done
for g in Cell Edge Face Outerface Side; do
    for t in Complex Double Float Integer; do
	${MT} $t.filenames ./tmp tbox tbox::Pointer pdat::${g}${t}ConstantRefine NDIM
    done
done

#
# basic patch data types and associated templates
#
for t in bool char dcomplex double int float; do
    ${MT} $t.filenames ./tmp pdat pdat::ArrayData NDIM\,$t
    ${MT} $t.filenames ./tmp pdat pdat::CellData NDIM\,$t
    ${MT} $t.filenames ./tmp tbox tbox::Array tbox::Pointer pdat::ArrayData NDIM\,$t
    ${MT} $t.filenames ./tmp tbox tbox::Pointer pdat::ArrayData  NDIM\,$t
    for g in Cell Edge Face Node Outeredge Outernode Outerface Outerside Side; do
	${MT} $t.filenames ./tmp pdat pdat::${g}Data NDIM\,$t
	${MT} $t.filenames ./tmp pdat pdat::${g}DataFactory NDIM\,$t
	${MT} $t.filenames ./tmp pdat pdat::${g}Variable NDIM\,$t
	${MT} $t.filenames ./tmp tbox tbox::Pointer pdat::${g}Data NDIM\,$t
	${MT} $t.filenames ./tmp tbox tbox::Pointer pdat::${g}DataFactory NDIM\,$t
	${MT} $t.filenames ./tmp tbox tbox::Pointer pdat::${g}Variable NDIM\,$t
    done
done

#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticXd ./tmp/*.C

sh ../../scripts/object.sh ./tmp automaticXd
sh ../../scripts/depend

rm -rf ./tmp


exit 0





