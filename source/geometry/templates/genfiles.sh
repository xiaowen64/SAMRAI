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
SkeletonPatchGeometry SkeletonGridGeometry SkeletonCoarsen \
SkeletonRefine CartesianPatchGeometry CartesianGridGeometry \
CartesianCellFloatLinearRefine CartesianCellComplexLinearRefine \
CartesianCellDoubleLinearRefine CartesianCellComplexWeightedAverage \
CartesianCellFloatWeightedAverage \
CartesianCellFloatConservativeLinearRefine \
CartesianCellComplexConservativeLinearRefine \
CartesianCellDoubleConservativeLinearRefine \
CartesianCellDoubleWeightedAverage \
CartesianOutersideDoubleWeightedAverage \
CartesianSideFloatConservativeLinearRefine \
CartesianSideComplexWeightedAverage CartesianSideDoubleWeightedAverage \
CartesianSideDoubleConservativeLinearRefine \
CartesianSideFloatWeightedAverage \
CartesianEdgeFloatConservativeLinearRefine \
CartesianEdgeComplexWeightedAverage CartesianEdgeDoubleWeightedAverage \
CartesianEdgeDoubleConservativeLinearRefine \
CartesianEdgeFloatWeightedAverage CartesianNodeFloatLinearRefine \
CartesianNodeComplexLinearRefine CartesianNodeDoubleLinearRefine \
CartesianOuterfaceFloatWeightedAverage \
CartesianOuterfaceComplexWeightedAverage \
CartesianOuterfaceDoubleWeightedAverage \
CartesianFaceComplexWeightedAverage CartesianFaceFloatWeightedAverage \
CartesianFaceFloatConservativeLinearRefine \
CartesianFaceDoubleWeightedAverage \
CartesianFaceDoubleConservativeLinearRefine 
do
${MT} default.filenames ./tmp geom $t NDIM
done

${MT} default.filenames ./tmp tbox tbox::Pointer geom::CartesianGridGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer geom::CartesianPatchGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer geom::SkeletonGridGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer geom::SkeletonPatchGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer geom::SkeletonGridGeometry NDIM

for g in Cell; do
    for t in Complex Double Float; do
	${MT} $t.filenames ./tmp tbox tbox::Pointer geom::Cartesian${g}${t}ConservativeLinearRefine NDIM
    done
done
for g in Edge Side; do
    for t in Double Float; do
	${MT} $t.filenames ./tmp tbox tbox::Pointer geom::Cartesian${g}${t}ConservativeLinearRefine NDIM
    done
done
for g in Cell Node; do
    for t in Complex Double Float; do
	${MT} $t.filenames ./tmp tbox tbox::Pointer geom::Cartesian${g}${t}LinearRefine NDIM
    done
done
for g in Cell Edge Face Outerface Side; do
    for t in Complex Double Float; do
	${MT} $t.filenames ./tmp tbox tbox::Pointer geom::Cartesian${g}${t}WeightedAverage NDIM
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





