//
// File:        $URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/solvers/vectors/SAMRAIVectorReal.I $
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 1704 $
// Modified:    $LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description: Vector class for real data on SAMRAI hierarchy.
//


#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace solv {

template<int DIM, class TYPE> inline
void SAMRAIVectorReal<DIM,TYPE>::setOutputStream(std::ostream& s)
{
   d_output_stream = &s;
}

template<int DIM, class TYPE> inline
std::ostream& SAMRAIVectorReal<DIM,TYPE>::getOutputStream()
{
   return(*d_output_stream);
}

template<int DIM, class TYPE> inline
const std::string& SAMRAIVectorReal<DIM,TYPE>::getName() const
{
   return(d_vector_name);
}

template<int DIM, class TYPE> inline
tbox::Pointer< hier::PatchHierarchy<DIM> > 
SAMRAIVectorReal<DIM,TYPE>::getPatchHierarchy() const
{
   return(d_hierarchy);
}

template<int DIM, class TYPE> inline
int SAMRAIVectorReal<DIM,TYPE>::getCoarsestLevelNumber() const
{
   return(d_coarsest_level);
}

template<int DIM, class TYPE> inline
int SAMRAIVectorReal<DIM,TYPE>::getFinestLevelNumber() const
{
   return(d_finest_level);
}

template<int DIM, class TYPE> inline
int SAMRAIVectorReal<DIM,TYPE>::getNumberOfComponents() const
{
   return(d_number_components);
}

template<int DIM, class TYPE> inline 
tbox::Pointer< hier::PatchData<DIM> > 
SAMRAIVectorReal<DIM,TYPE>::getComponentPatchData(
   const int comp_id, 
   const hier::Patch<DIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components); 
#endif
   return( patch.getPatchData(d_component_data_id[comp_id]) ); 
}

template<int DIM, class TYPE> inline 
tbox::Pointer< hier::PatchData<DIM> > 
SAMRAIVectorReal<DIM,TYPE>::getComponentPatchData(
   const tbox::Pointer< hier::Variable<DIM> >& var,
   const hier::Patch<DIM>& patch) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(!var.isNull());
   TBOX_ASSERT(d_variableid_2_vectorcomponent_map[
          var->getInstanceIdentifier()] >= 0);
#endif
   return( patch.getPatchData(
                 d_component_data_id[
                    d_variableid_2_vectorcomponent_map[
                     var->getInstanceIdentifier()] ]) );
}

template<int DIM, class TYPE> inline
tbox::Pointer< hier::Variable<DIM> >
SAMRAIVectorReal<DIM,TYPE>::getComponentVariable(
   const int comp_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
#endif
   return(d_component_variable[comp_id]);
   
}

template<int DIM, class TYPE> inline 
int SAMRAIVectorReal<DIM,TYPE>::getComponentDescriptorIndex(
   const int comp_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
#endif
   return( d_component_data_id[comp_id] ); 
}

template<int DIM, class TYPE> inline
int SAMRAIVectorReal<DIM,TYPE>::getControlVolumeIndex(
   const int comp_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   TBOX_ASSERT(comp_id >= 0 && comp_id < d_number_components);
#endif
   return( d_control_volume_data_id[comp_id] );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
