##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/patchdata/templates/genfiles.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 2153 $
## Modified:    $LastChangedDate: 2008-04-25 14:08:49 -0700 (Fri, 25 Apr 2008) $
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
NodeIndex NodeIterator \
FaceIndex FaceIterator \
OuternodeDoubleConstantCoarsen;
do
  ${MT} default.filenames ./tmpXd pdat $t NDIM
done

for g in Cell Side Outerside Edge Node Outerface Face; do
    for t in Complex Double Float; do
	${MT} ${t}.filenames ./tmpXd pdat ${g}${t}LinearTimeInterpolateOp NDIM
    done
done

for g in Cell Side Edge Outerface Face; do 
    for t in Complex Double Float Integer; do
	${MT} ${t}.filenames ./tmpXd pdat ${g}${t}ConstantRefine NDIM
    done
done

for g in Node; do 
    for t in Complex Double Float Integer; do
	${MT} ${t}.filenames ./tmpXd pdat ${g}${t}Injection NDIM
    done
done


# NonDim
for t in bool char dcomplex double int float; do
   ${MT} $t.filenames ./tmpNond pdat pdat::CopyOperation $t
done
for t in dcomplex double int float; do
   ${MT} $t.filenames ./tmpNond pdat pdat::SumOperation $t
done


for g in Node; do
    for t in Complex Double Float Integer; do
	${MT} $t.filenames ./tmpXd tbox tbox::Pointer pdat::${g}${t}Injection NDIM
    done
done
for g in Outernode; do
    for t in Double; do
	${MT} $t.filenames ./tmpXd tbox tbox::Pointer pdat::${g}${t}ConstantCoarsen NDIM
    done
done
for g in Cell Edge Face Outerface Side; do
    for t in Complex Double Float Integer; do
	${MT} $t.filenames ./tmpXd tbox tbox::Pointer pdat::${g}${t}ConstantRefine NDIM
    done
done

for t in double float bool dcomplex int; do
    ${MT} $t.filenames ./tmpXd pdat pdat::MBDataUtilities NDIM\,$t
    ${MT} $t.filenames ./tmpXd pdat pdat::MultiblockCellDataTranslator NDIM\,$t
    ${MT} $t.filenames ./tmpXd pdat pdat::MultiblockEdgeDataTranslator NDIM\,$t
    ${MT} $t.filenames ./tmpXd pdat pdat::MultiblockFaceDataTranslator NDIM\,$t
    ${MT} $t.filenames ./tmpXd pdat pdat::MultiblockNodeDataTranslator NDIM\,$t
    ${MT} $t.filenames ./tmpXd pdat pdat::MultiblockSideDataTranslator NDIM\,$t
done

#
# basic patch data types and associated templates
#
for t in bool char dcomplex double int float; do
    ${MT} $t.filenames ./tmpXd pdat pdat::ArrayData NDIM\,$t
    ${MT} $t.filenames ./tmpXd pdat pdat::CellData NDIM\,$t
    ${MT} $t.filenames ./tmpXd tbox tbox::Array tbox::Pointer pdat::ArrayData NDIM\,$t
    ${MT} $t.filenames ./tmpXd tbox tbox::Pointer pdat::ArrayData  NDIM\,$t
    for g in Cell Edge Face Node Outeredge Outernode Outerface Outerside Side; do
	${MT} $t.filenames ./tmpXd pdat pdat::${g}Data NDIM\,$t
	${MT} $t.filenames ./tmpXd pdat pdat::${g}DataFactory NDIM\,$t
	${MT} $t.filenames ./tmpXd pdat pdat::${g}Variable NDIM\,$t
	${MT} $t.filenames ./tmpXd tbox tbox::Pointer pdat::${g}Data NDIM\,$t
	${MT} $t.filenames ./tmpXd tbox tbox::Pointer pdat::${g}DataFactory NDIM\,$t
	${MT} $t.filenames ./tmpXd tbox tbox::Pointer pdat::${g}Variable NDIM\,$t
    done
done
for t in bool char dcomplex double int float; do
    for op in CopyOperation; do
       ${MT} $t.filenames ./tmpXd pdat pdat::ArrayDataOperationUtilities NDIM\,$t\,${op} $t
    done
done
for t in bool char dcomplex double int float; do
    for op in SumOperation; do
       ${MT} $t.filenames ./tmpXd pdat pdat::ArrayDataOperationUtilities NDIM\,$t\,${op} $t
    done
done


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





