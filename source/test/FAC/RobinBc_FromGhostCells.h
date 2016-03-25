/*
  File:		$RCSfile$
  Copyright:	(c) 1997-2005 The Regents of the University of California
  Revision:	$Revision: 173 $
  Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
  Description:	RobinBc_FromGhostCells class declaration
*/

#ifndef included_RobinBc_FromGhostCells
#define included_RobinBc_FromGhostCells


#ifndef included_tbox_ArrayData
#include "ArrayData.h"
#endif

#ifndef included_solv::CartesianRobinBc
#include "CartesianRobinBc.h"
#endif

#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif

#ifndef included_String
#include "tbox/String.h"
#endif


using namespace SAMRAI;


/*!
  @brief Implementation of solv::CartesianRobinBc,
  specializing in setting bc to given ghost cell values.

  This implementation of solv::CartesianRobinBc is a handy
  class for specifying Robin boundary conditions when you
  have ghost cell data stored.  It sets the coefficients
  a and g (see solv::CartesianRobinBc) in such a way that
  the Robin boundary condition specifies that the value
  at the ghost cell center is the value you stored.

  This class specializes to cell data, hence it does not
  override any solv::CartesianRobinBc functions that does
  not apply to cell data.
*/
template<int NDIM>
class RobinBc_FromGhostCells :
  public solv::CartesianRobinBc<NDIM>
{

public:

  /*!
    @brief Constructor.

    If you want standard output and logging,
    pass in valid pointers for those streams.

    @param object_name Object name
    @param out_stream tbox::Pointer to output stream (NULL if unused).
    @param log_stream tbox::Pointer to .og stream (NULL if unused).
  */
  RobinBc_FromGhostCells(
    const string &object_name ,
    ostream *out_stream=NULL ,
    ostream *log_stream=NULL );



  //@{

  /*!
    @name Functions supporting setBcCoefOnBoundaryBoxAtSides
  */

  void setBcCoefOnBoundaryBoxAtSides (
    tbox::Pointer<pdat::ArrayData<NDIM,double> > &acoef_data ,
    tbox::Pointer<pdat::ArrayData<NDIM,double> > &gcoef_data ,
    const hier::Patch<NDIM> &patch ,
    const hier::BoundaryBox<NDIM> &bdry_box ,
    double fill_time ) const;

  //@}

   /*!
     @brief Set descriptor index to use to get ghost cell data.

     The index @em must correspond to cell-centered double data.

     @internal If need be, the descriptor can be allowed to point
     to outerside data.  This may be useful if the user must work
     with data that does not have ghost cells.  Allocating
     outerside data is less expensive than allocating cell data
     with ghost cells.  This feature is implemented but not yet
     tested, so it is not officially supported.
   */
   void setGhostDataId( int ghost_data_id );

  //@{
private:
  /*!
    @name Private state variables for solution.
  */

  /*!
    @brief hier::Index of variable to use for ghost cell data.
  */
  int d_ghost_data_id;

  //@}



  //@{
private:

   /*!
     @brief Object name.
   */
   string d_object_name;
  /*!
    @name Output streams.
  */
  /*!
    @brief Output stream pointer.

    If set to NULL, no output.
  */
  ostream *d_ostream;

  /*!
    @brief Log stream pointer.

    If set to NULL, no logging.
  */
  ostream *d_lstream;
  //@}


};


#endif	// included_RobinBc_FromGhostCells
