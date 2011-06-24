/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Vector class for real data on SAMRAI hierarchy. 
 *
 ************************************************************************/

namespace SAMRAI {
namespace solv {

template<class TYPE>
SAMRAI_INLINE_KEYWORD
void SAMRAIVectorReal<TYPE>::setOutputStream(
   std::ostream& s)
{
   d_output_stream = &s;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
std::ostream& SAMRAIVectorReal<TYPE>::getOutputStream()
{
   return *d_output_stream;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
const std::string& SAMRAIVectorReal<TYPE>::getName() const
{
   return d_vector_name;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
tbox::Pointer<hier::PatchHierarchy>
SAMRAIVectorReal<TYPE>::getPatchHierarchy() const
{
   return d_hierarchy;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int SAMRAIVectorReal<TYPE>::getCoarsestLevelNumber() const
{
   return d_coarsest_level;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int SAMRAIVectorReal<TYPE>::getFinestLevelNumber() const
{
   return d_finest_level;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int SAMRAIVectorReal<TYPE>::getNumberOfComponents() const
{
   return d_number_components;
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
tbox::Pointer<hier::PatchData>
SAMRAIVectorReal<TYPE>::getComponentPatchData(
   const int comp_id,
   const hier::Patch& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
#endif
   return patch.getPatchData(d_component_data_id[comp_id]);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
tbox::Pointer<hier::PatchData>
SAMRAIVectorReal<TYPE>::getComponentPatchData(
   const tbox::Pointer<hier::Variable>& var,
   const hier::Patch& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!var.isNull());
   TBOX_ASSERT(d_variableid_2_vectorcomponent_map[
         var->getInstanceIdentifier()] >= 0);
#endif
   return patch.getPatchData(
      d_component_data_id[
         d_variableid_2_vectorcomponent_map[
            var->getInstanceIdentifier()]]);
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
tbox::Pointer<hier::Variable>
SAMRAIVectorReal<TYPE>::getComponentVariable(
   const int comp_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
#endif
   return d_component_variable[comp_id];

}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int SAMRAIVectorReal<TYPE>::getComponentDescriptorIndex(
   const int comp_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
#endif
   return d_component_data_id[comp_id];
}

template<class TYPE>
SAMRAI_INLINE_KEYWORD
int SAMRAIVectorReal<TYPE>::getControlVolumeIndex(
   const int comp_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
#endif
   return d_control_volume_data_id[comp_id];
}

}
}
