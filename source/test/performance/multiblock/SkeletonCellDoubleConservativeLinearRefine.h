/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Conservative linear refine operator for cell-centered
 *                double data on a Skeleton mesh. 
 *
 ************************************************************************/

#ifndef included_SkeletonCellDoubleConservativeLinearRefineXD
#define included_SkeletonCellDoubleConservativeLinearRefineXD

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/tbox/Pointer.h"
#ifndef included_String
#include <string>
#define included_String
#endif
#include "SAMRAI/hier/RefineOperator.h"

using namespace std;
using namespace SAMRAI;

/**
 * Class SkeletonCellDoubleConservativeLinearRefine implements
 * conservative linear interpolation for cell-centered double patch data
 * defined over a Skeleton mesh.  It is derived from the base class
 * hier::RefineOperator.  The numerical operations for the interpolation
 * use FORTRAN numerical routines.
 *
 * The findRefineOperator() operator function returns true if the input
 * variable is cell-centered double, and the string is
 * "CONSERVATIVE_LINEAR_REFINE".
 *
 * @see hier::RefineOperator
 */

class SkeletonCellDoubleConservativeLinearRefine:
   public hier::RefineOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   SkeletonCellDoubleConservativeLinearRefine(
      const tbox::Dimension& dim);

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SkeletonCellDoubleConservativeLinearRefine();

   /**
    * Return true if the variable and name string match cell-centered
    * double conservative linear interpolation; otherwise, return false.
    */
   bool
   findRefineOperator(
      const tbox::Pointer<hier::Variable>& var,
      const string& op_name) const;

   /**
    * The priority of cell-centered double conservative linear is 0.
    * It will be performed before any user-defined interpolation operations.
    */
   int
   getOperatorPriority() const;

   /**
    * The stencil width of the conservative linear interpolation operator is
    * the vector of ones.
    */
   hier::IntVector
   getStencilWidth() const;

   /**
    * Refine the source component on the coarse patch to the destination
    * component on the fine patch using the cell-centered double conservative
    * linear interpolation operator.  Interpolation is performed on the
    * intersection of the destination patch and the fine box.   It is assumed
    * that the coarse patch contains sufficient data for the stencil width of
    * the refinement operator.
    */
   void
   refine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const int dst_component,
      const int src_component,
      const hier::BoxOverlap& fine_overlap,
      const hier::IntVector& ratio) const;

   void
   refine(
      hier::Patch& fine,
      const hier::Patch& coarse,
      const int dst_component,
      const int src_component,
      const hier::Box& fine_box,
      const hier::IntVector& ratio) const;

   /**
    * Set the dx, the distance between mesh nodes.
    */
   void
   setDx(
      const int level_number,
      const double* dx);

private:
   /**
    * Return the dx
    */
   void
   getDx(
      const int level_number,
      double* dx) const;

   tbox::Array<tbox::Array<double> > d_dx;

};

#endif
