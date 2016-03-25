//
// File:	SkeletonCellDoubleWeightedAverage.h
// Package:	SAMRAI geometry
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description: Weighted averaging operator for cell-centered double data on 
//              a Skeleton mesh.
//

#ifndef included_SkeletonCellDoubleWeightedAverageXD
#define included_SkeletonCellDoubleWeightedAverageXD

#include "SAMRAI_config.h"
#include "Box.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "CoarsenOperator.h"

using namespace std;
using namespace SAMRAI;

/**
 * Class SkeletonCellDoubleWeightedAverageX implements conservative
 * cell-weighted averaging for cell-centered double patch data defined over a 
 * Skeleton mesh.  It is derived from the xfer::CoarsenOperator<NDIM> base class.
 * The numerical operations for the averaging use FORTRAN numerical routines.
 *
 * The findCoarsenOperator() operator function returns true if the input 
 * variable is cell-centered double, and the string is "CONSERVATIVE_COARSEN".
 * 
 * @see xfer::CoarsenOperator
 */

class SkeletonCellDoubleWeightedAverage 
: public xfer::CoarsenOperator<NDIM>
{
public:
   /**
    * Uninteresting default constructor.
    */
   SkeletonCellDoubleWeightedAverage();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SkeletonCellDoubleWeightedAverage();

   /**
    * Return true if the variable and name string match cell-centered 
    * double weighted averaging; otherwise, return false.
    */
   bool findCoarsenOperator(const tbox::Pointer< hier::Variable<NDIM> >& var,
                            const string &op_name) const; 

   /**
    * Return name string identifier of this coarsening operation.
    */
   const string& getOperatorName() const;

   /**
    * The priority of cell-centered double weighted averaging is 0.
    * It will be performed before any user-defined coarsen operations. 
    */
   int getOperatorPriority() const;

   /**
    * The stencil width of the weighted averaging operator is the vector of 
    * zeros.  That is, its stencil does not extend outside the fine box.
    */
   hier::IntVector<NDIM> getStencilWidth() const;

   /**
    * Coarsen the source component on the fine patch to the destination
    * component on the coarse patch using the cell-centered double weighted 
    * averaging operator.  Coarsening is performed on the intersection of 
    * the destination patch and the coarse box.  It is assumed that the 
    * fine patch contains sufficient data for the stencil width of the 
    * coarsening operator.
    */
   void coarsen(hier::Patch<NDIM>& coarse,
                const hier::Patch<NDIM>& fine,
                const int dst_component,
                const int src_component,
                const hier::Box<NDIM>& coarse_box,
                const hier::IntVector<NDIM>& ratio) const;

   /**
    * Set the dx, the distance between mesh nodes.  
    */
   void setDx(const int level_number,
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

