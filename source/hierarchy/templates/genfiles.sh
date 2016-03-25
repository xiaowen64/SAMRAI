##
## File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/hierarchy/templates/genfiles.sh $
## Package:     SAMRAI templates
## Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
## Revision:    $LastChangedRevision: 2040 $
## Modified:    $LastChangedDate: 2008-03-11 15:05:44 -0700 (Tue, 11 Mar 2008) $
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

for t in BasePatchHierarchy BasePatchLevel BinaryTree BoundaryBox \
    BoundaryBoxUtils \
    BoundaryLookupTable Box BoxArray BoxComm BoxGeometry BoxGraph \
    BoxGraphUtilities BoxIOUtility BoxList BoxOverlap BoxTop BoxTree \
    BoxTreeNode BoxUtilities CoarseFineBoundary GridGeometry Index IntVector \
    LayerEdgeSet LayerHierarchy LayerNode LayerNodeSet \
    LocallyActiveDataPatchLevelManager LocallyActiveVariableDatabase \
    MBUtilities MultiblockDataTranslator MultiblockGridGeometry \
    MultiblockPatchHierarchy \
    MultiblockPatchLevel Patch PatchConfigurationUtilities PatchData \
    PatchDataFactory PatchDescriptor PatchFactory PatchGeometry \
    PatchHierarchy PatchLevel PatchLevelFactory Variable VariableDatabase

do
  ${MT} default.filenames ./tmp hier $t NDIM
done

${MT} default.filenames ./tmp tbox tbox::Array hier::ComponentSelector

${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::ComponentSelector
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::VariableContext

${MT} default.filenames ./tmp tbox tbox::Array hier::BinaryTree\<NDIM\>::binTreeNode


${MT} default.filenames ./tmp tbox tbox::Array hier::BoundaryBox NDIM
${MT} default.filenames ./tmp tbox tbox::Array hier::Box NDIM
${MT} default.filenames ./tmp tbox tbox::Array hier::BoxArray NDIM
${MT} default.filenames ./tmp tbox tbox::Array hier::BoxGraph\<NDIM\>::GraphNode
${MT} default.filenames ./tmp tbox tbox::Array hier::BoxList NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::BoxTree NDIM
${MT} default.filenames ./tmp tbox tbox::Array hier::BoxTreeNode\<NDIM\>::Triple
${MT} default.filenames ./tmp tbox tbox::Array hier::CoarseFineBoundary NDIM
${MT} default.filenames ./tmp tbox tbox::Array hier::Index NDIM
${MT} default.filenames ./tmp tbox tbox::Array hier::IntVector NDIM
${MT} default.filenames ./tmp tbox tbox::Array hier::LayerEdgeSet NDIM
${MT} default.filenames ./tmp tbox tbox::Array hier::LayerNodeSet NDIM

${MT} default.filenames ./tmp tbox tbox::Array hier::PatchConfigurationUtilities\<NDIM\>::PatchLevelInfo
${MT} default.filenames ./tmp tbox tbox::Array hier::PatchConfigurationUtilities\<NDIM\>::PatchInfo*
${MT} default.filenames ./tmp tbox tbox::Array hier::PatchConfigurationUtilities\<NDIM\>::NeighborPatchInfo

${MT} default.filenames ./tmp tbox tbox::Array tbox::Array hier::Box NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array hier::BoxArray NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array hier::Index NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array hier::IntVector NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array hier::BoxList NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array hier::BoundaryBox NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Array hier::BoundaryBox NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Array hier::CoarseFineBoundary NDIM

${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Array tbox::Pointer hier::PatchLevel NDIM

${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Pointer hier::PatchLevel NDIM

${MT} default.filenames ./tmp tbox tbox::Array tbox::List hier::IntVector NDIM

${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::BoxOverlap NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::Index NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::GridGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::Patch NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::PatchData NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::PatchDataFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::PatchHierarchy NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::PatchLevel NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::Variable NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::LocallyActiveDataPatchLevelManager NDIM

${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::PatchConfigurationUtilities NDIM

${MT} default.filenames ./tmp tbox tbox::List hier::BoundaryBox NDIM
${MT} default.filenames ./tmp tbox tbox::List hier::BoxTreeNode\<NDIM\>::Triple
${MT} default.filenames ./tmp tbox tbox::List hier::Box NDIM
${MT} default.filenames ./tmp tbox tbox::List hier::Index NDIM
${MT} default.filenames ./tmp tbox tbox::List hier::IntVector NDIM

${MT} default.filenames ./tmp tbox tbox::Pointer hier::ComponentSelector
${MT} default.filenames ./tmp tbox tbox::Pointer hier::VariableContext

${MT} default.filenames ./tmp tbox tbox::Pointer hier::BasePatchHierarchy NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BasePatchLevel NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BinaryTree NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BoxGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BoxGraph NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BoxList NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BoxOverlap NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BoxTop NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BoxTree NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::BoxTreeNode NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::CoarseFineBoundary NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::GridGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::Index NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::Patch NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchData NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchDataFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchDescriptor NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchHierarchy NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchLevel NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchLevelFactory NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::Variable NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::LocallyActiveDataPatchLevelManager NDIM

${MT} default.filenames ./tmp tbox tbox::Pointer hier::PatchConfigurationUtilities NDIM

${MT} default.filenames ./tmp tbox tbox::Array tbox::Array tbox::Pointer hier::MultiblockPatchLevel NDIM
${MT} default.filenames ./tmp tbox tbox::Array tbox::List hier::MultiblockPatchHierarchy\<NDIM\>::Neighbor
${MT} default.filenames ./tmp tbox tbox::Array tbox::Pointer hier::MultiblockPatchLevel NDIM

${MT} default.filenames ./tmp tbox tbox::List hier::MultiblockPatchHierarchy\<NDIM\>::Neighbor

${MT} default.filenames ./tmp tbox tbox::Pointer hier::MultiblockGridGeometry NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::MultiblockPatchHierarchy NDIM
${MT} default.filenames ./tmp tbox tbox::Pointer hier::MultiblockPatchLevel NDIM

#
# other list templates
#

${MT} default.filenames ./tmp tbox tbox::List tbox::Pointer hier::Variable NDIM

#
# now copy the new template files into the repository
#

sh ../../scripts/copy-if-change ./automaticXd ./tmp/*.C

sh ../../scripts/object.sh ./tmp automaticXd
sh ../../scripts/depend

rm -rf ./tmp


exit 0

