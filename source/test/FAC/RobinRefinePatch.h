/*
  File:		$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/test/FAC/RobinRefinePatch.h $
  Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
  Revision:	$LastChangedRevision: 1704 $
  Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
  Description:	HyprePoisson class declaration
*/

#ifndef included_RobinRefinePatch
#define included_RobinRefinePatch


#include <string>

#include "tbox/Pointer.h"
#include "tbox/Database.h"


#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif


/*
  SAMRAI classes
*/

#ifndef solv::CartesianRobinBcHelper<NDIM>
#include "CartesianRobinBcHelper.h"
#endif

#ifndef solv::RobinBcCoefStrategy<NDIM>
#include "RobinBcCoefStrategy.h"
#endif

#ifndef xfer::RefinePatchStrategy<NDIM>
#include "RefinePatchStrategy.h"
#endif


#ifndef LACKS_NAMESPACE
namespace SAMRAI {
#endif


/*!
  @brief Implements xfer::RefinePatchStrategy<NDIM> for specific case of
  using Robin boundary conditions on cell-centered data.

  This class inherits and implements virtual functions from
  - xfer::RefinePatchStrategy<NDIM>

  It uses the helper class solv::CartesianRobinBcHelper<NDIM> to help
  fill the ghost cells.  It requires an implementation of
  solv::RobinBcCoefStrategy<NDIM> to give it the boundary condition
  coefficients.
*/
template<int NDIM>
class RobinRefinePatch :
  public xfer::RefinePatchStrategy<NDIM>
{

public:

  /*!
    @brief Constructor.
    @param object_name Ojbect name
    @param database Input database
    @param out_stream Standard output stream
    @param log_stream Log output stream
  */
  RobinRefinePatch( const string &object_name=string(),
		     solv::RobinBcCoefStrategy<NDIM> *coef_strategy=NULL );



   /*!
     @brief Destructor.
   */
   ~RobinRefinePatch();



   //@{ @name xfer::RefinePatchStrategy virtuals

   virtual void setPhysicalBoundaryConditions (
      hier::Patch<NDIM> &patch ,
      const double fill_time ,
      const hier::IntVector<NDIM> &ghost_width_to_fill );
   hier::IntVector<NDIM> getRefineOpStencilWidth () const;
   virtual void preprocessRefineBoxes (
      hier::Patch<NDIM> &fine ,
      const hier::Patch<NDIM> &coarse ,
      const hier::BoxList<NDIM> &fine_boxes ,
      const hier::IntVector<NDIM> &ratio );
   virtual void preprocessRefine (
      hier::Patch<NDIM> &fine ,
      const hier::Patch<NDIM> &coarse ,
      const hier::Box<NDIM> &fine_box ,
      const hier::IntVector<NDIM> &ratio );
   virtual void postprocessRefineBoxes (
      hier::Patch<NDIM> &fine ,
      const hier::Patch<NDIM> &coarse ,
      const hier::BoxList<NDIM> &fine_box ,
      const hier::IntVector<NDIM> &ratio );
   virtual void postprocessRefine (
      hier::Patch<NDIM> &fine ,
      const hier::Patch<NDIM> &coarse ,
      const hier::Box<NDIM> &fine_boxes ,
      const hier::IntVector<NDIM> &ratio );

   //@}



   /*!
     @brief Provide an implementation of the solv::RobinBcCoefStrategy
     for determining the boundary coefficients.

     Provide the implementation that can be used to set the
     Robin bc coefficients.

     Use of this function excludes the use of setCoefForParallelpipedDomain().

     @param coef_strategy tbox::Pointer to a concrete inmplementation of
            the coefficient strategy.
   */
   void setCoefImplementation(
      const solv::RobinBcCoefStrategy<NDIM> *coef_strategy );


   /*!
     @brief Set the data id that should be filled when setting
     physical boundary conditions.
   */
   void setDataId( int data_id );



   /*!
     @brief Set whether boundary filling should assume homogeneous
     conditions.
   */
   void setHomogeneousBc( bool is_homogeneous );



private:


   string d_name;

   solv::CartesianRobinBcHelper<NDIM> d_bc_filler;

   int d_data_id;

   bool d_homogeneous_bc;


};


#ifndef LACKS_NAMESPACE
}
#endif

#ifndef DEBUG_NO_INLINE
// #include "RobinRefinePatch.I"
#endif

#endif	// included_HyprePoisson
