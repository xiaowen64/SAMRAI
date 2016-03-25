##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/geometry/templates/genfiles.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 2158 $
## Modified:    $LastChangedDate: 2008-04-28 17:04:43 -0700 (Mon, 28 Apr 2008) $
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
BlockPatchGeometry BlockGridGeometry SkeletonCoarsen \
SkeletonPatchGeometry SkeletonGridGeometry SkeletonCoarsen \
SkeletonRefine CartesianPatchGeometry CartesianGridGeometry \
CartesianCellDoubleLinearRefine \
CartesianCellDoubleConservativeLinearRefine \
CartesianCellDoubleWeightedAverage \
CartesianOutersideDoubleWeightedAverage \
CartesianSideDoubleWeightedAverage \
CartesianSideDoubleConservativeLinearRefine \
CartesianEdgeDoubleWeightedAverage \
CartesianEdgeDoubleConservativeLinearRefine \
CartesianNodeDoubleLinearRefine \
CartesianOuterfaceDoubleWeightedAverage \
CartesianFaceDoubleWeightedAverage \
CartesianFaceDoubleConservativeLinearRefine 
do
${MT} default.filenames ./tmp geom $t NDIM
done

for t in \
CartesianCellFloatLinearRefine \
CartesianCellFloatWeightedAverage \
CartesianCellFloatConservativeLinearRefine \
CartesianSideFloatConservativeLinearRefine \
CartesianSideFloatWeightedAverage \
CartesianEdgeFloatConservativeLinearRefine \
CartesianEdgeFloatWeightedAverage \
CartesianNodeFloatLinearRefine \
CartesianOuterfaceFloatWeightedAverage \
CartesianFaceFloatWeightedAverage \
CartesianFaceFloatConservativeLinearRefine 
do
${MT} Float.filenames ./tmp geom $t NDIM
done

for t in \
CartesianCellComplexLinearRefine \
CartesianCellComplexWeightedAverage \
CartesianCellComplexConservativeLinearRefine \
CartesianSideComplexWeightedAverage \
CartesianEdgeComplexWeightedAverage \
CartesianNodeComplexLinearRefine \
CartesianOuterfaceComplexWeightedAverage \
CartesianFaceComplexWeightedAverage CartesianFaceFloatWeightedAverage
do
${MT} Complex.filenames ./tmp geom $t NDIM
done



${MT} default.filenames ./tmp tbox tbox::Pointer geom::BlockGridGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer geom::BlockPatchGeometry NDIM
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





