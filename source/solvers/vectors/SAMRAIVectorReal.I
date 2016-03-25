//
// File:        SAMRAIVectorReal.I
// Package:     SAMRAI solvers
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Vector class for real data on SAMRAI hierarchy.
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_assert
#define included_assert
#include <assert.h>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace SAMRAI {
    namespace solv {

template<int DIM, class TYPE> inline
void SAMRAIVectorReal<DIM,TYPE>::setOutputStream(ostream& s)
{
   d_output_stream = &s;
}

template<int DIM, class TYPE> inline
ostream& SAMRAIVectorReal<DIM,TYPE>::getOutputStream()
{
   return(*d_output_stream);
}

template<int DIM, class TYPE> inline
const string& SAMRAIVectorReal<DIM,TYPE>::getName() const
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
   assert(comp_id >= 0 && comp_id < d_number_components); 
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
   assert(!var.isNull());
   assert(d_variableid_2_vectorcomponent_map[
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
   assert(comp_id >= 0 && comp_id < d_number_components);
#endif
   return(d_component_variable[comp_id]);
   
}

template<int DIM, class TYPE> inline 
int SAMRAIVectorReal<DIM,TYPE>::getComponentDescriptorIndex(
   const int comp_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(comp_id >= 0 && comp_id < d_number_components);
#endif
   return( d_component_data_id[comp_id] ); 
}

template<int DIM, class TYPE> inline
int SAMRAIVectorReal<DIM,TYPE>::getControlVolumeIndex(
   const int comp_id) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(comp_id >= 0 && comp_id < d_number_components);
#endif
   return( d_control_volume_data_id[comp_id] );
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}
