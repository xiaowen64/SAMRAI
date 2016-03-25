/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/RobinRefinePatch.C $
  Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	HyprePoisson class implementation
*/

#include "SAMRAI_config.h"

#include "RobinRefinePatch.h"

#include IOMANIP_HEADER_FILE


#ifdef DEBUG_NO_INLINE
// #include "RobinRefinePatch.I"
#endif

#ifndef LACKS_NAMESPACE
namespace SAMRAI {
#endif

template<int NDIM>
RobinRefinePatch<NDIM>::RobinRefinePatch(
   const string &object_name,
   solv::RobinBcCoefStrategy<NDIM> *coef_strategy)
   :
   d_name(object_name),
   d_bc_filler(),
   d_data_id(-1),
   d_homogeneous_bc(false)
{
   if ( coef_strategy ) {
      d_bc_filler.setCoefImplementation( coef_strategy );
   }
   return;
}

template<int NDIM>
RobinRefinePatch<NDIM>::~RobinRefinePatch() {
   return;
}

template<int NDIM>
void RobinRefinePatch<NDIM>::setCoefImplementation (
   const solv::RobinBcCoefStrategy<NDIM> *coef_strategy )
{
   d_bc_filler.setCoefImplementation(coef_strategy);
   return;
}

template<int NDIM>
void RobinRefinePatch<NDIM>::setDataId( int data_id ) {
   d_data_id = data_id;
   return;
}

template<int NDIM>
void RobinRefinePatch<NDIM>::setHomogeneousBc( bool is_homogeneous ) {
   d_homogeneous_bc = is_homogeneous;
}

/*
***********************************************************************

 Virtual functions or xfer::RefinePatchStrategy<NDIM>.

************************************************************************
*/

template<int NDIM>
void RobinRefinePatch<NDIM>::setPhysicalBoundaryConditions (
   hier::Patch<NDIM> &patch ,
   const double fill_time ,
   const hier::IntVector<NDIM> &ghost_width_to_fill ) {

   d_bc_filler.setBoundaryValuesInCells( patch ,
					 0.0 ,
					 ghost_width_to_fill ,
					 d_data_id ,
					 d_homogeneous_bc );
   return;
}

template<int NDIM>
hier::IntVector<NDIM> RobinRefinePatch<NDIM>::getRefineOpStencilWidth() const {
   return hier::IntVector<NDIM>(1);
}

template<int NDIM>
void RobinRefinePatch<NDIM>::preprocessRefineBoxes (
      hier::Patch<NDIM> &fine ,
      const hier::Patch<NDIM> &coarse ,
      const hier::BoxList<NDIM> &fine_boxes ,
      const hier::IntVector<NDIM> &ratio ) {
   return;
}
template<int NDIM>
void RobinRefinePatch<NDIM>::preprocessRefine (
      hier::Patch<NDIM> &fine ,
      const hier::Patch<NDIM> &coarse ,
      const hier::Box<NDIM> &fine_box ,
      const hier::IntVector<NDIM> &ratio ) {
   return;
}
template<int NDIM>
void RobinRefinePatch<NDIM>::postprocessRefineBoxes (
      hier::Patch<NDIM> &fine ,
      const hier::Patch<NDIM> &coarse ,
      const hier::BoxList<NDIM> &fine_box ,
      const hier::IntVector<NDIM> &ratio ) {
   return;
}
template<int NDIM>
void RobinRefinePatch<NDIM>::postprocessRefine (
      hier::Patch<NDIM> &fine ,
      const hier::Patch<NDIM> &coarse ,
      const hier::Box<NDIM> &fine_boxes ,
      const hier::IntVector<NDIM> &ratio ) {
   return;
}

#ifndef LACKS_NAMESPACE
}
#endif

