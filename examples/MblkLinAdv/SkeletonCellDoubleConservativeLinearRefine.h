//
// File:	SkeletonCellDoubleConservativeLinearRefine.h
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Conservative linear refine operator for cell-centered 
//              double data on a Skeleton mesh.
//

#ifndef included_SkeletonCellDoubleConservativeLinearRefineXD
#define included_SkeletonCellDoubleConservativeLinearRefineXD


#include "SAMRAI_config.h"
#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "RefineOperator.h"

using namespace std;
using namespace SAMRAI;

/**
 * Class SkeletonCellDoubleConservativeLinearRefine implements
 * conservative linear interpolation for cell-centered double patch data 
 * defined over a Skeleton mesh.  It is derived from the base class 
 * xfer::RefineOperator<DIM>.  The numerical operations for the interpolation 
 * use FORTRAN numerical routines.
 *
 * The findRefineOperator() operator function returns true if the input 
 * variable is cell-centered double, and the string is 
 * "CONSERVATIVE_LINEAR_REFINE".
 * 
 * @see xfer::RefineOperator
 */

class SkeletonCellDoubleConservativeLinearRefine 
: public xfer::RefineOperator<NDIM>
{
public:
   /**
    * Uninteresting default constructor.
    */
   SkeletonCellDoubleConservativeLinearRefine();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SkeletonCellDoubleConservativeLinearRefine();

   /**
    * Return true if the variable and name string match cell-centered 
    * double conservative linear interpolation; otherwise, return false.
    */
   bool findRefineOperator(const tbox::Pointer< hier::Variable<NDIM> >& var,
                           const string &op_name) const; 

   /**
    * Return name string identifier of this refinement operator.
    */
   const string& getOperatorName() const;

   /**
    * The priority of cell-centered double conservative linear is 0.
    * It will be performed before any user-defined interpolation operations. 
    */
   int getOperatorPriority() const;

   /**
    * The stencil width of the conservative linear interpolation operator is 
    * the vector of ones.
    */
   hier::IntVector<NDIM> getStencilWidth() const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the cell-centered double conservative
    * linear interpolation operator.  Interpolation is performed on the 
    * intersection of the destination patch and the fine box.   It is assumed 
    * that the coarse patch contains sufficient data for the stencil width of 
    * the refinement operator.
    */
   void refine(hier::Patch<NDIM>& fine,
               const hier::Patch<NDIM>& coarse,
               const int dst_component,
               const int src_component,
               const hier::Box<NDIM>& fine_box,
               const hier::IntVector<NDIM>& ratio) const;

   /**
    * Set the dx, the distance between mesh nodes.  
    */
   void setDx( const int level_number,
               const double* dx);

private:

   /**
    * Return the dx  
    */
   void getDx( const int level_number,
               double* dx) const;

   string d_name_id;
   tbox::Array<tbox::Array<double> > d_dx;

};


#endif
