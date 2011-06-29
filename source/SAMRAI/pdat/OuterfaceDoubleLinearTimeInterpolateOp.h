/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for double outerface patch data. 
 *
 ************************************************************************/

#ifndef included_pdat_OuterfaceDoubleLinearTimeInterpolateOp
#define included_pdat_OuterfaceDoubleLinearTimeInterpolateOp

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/tbox/Pointer.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/**
 * Class OuterfaceDoubleLinearTimeInterpolateOp implements standard
 * linear time interpolation for double outerface patch data. Recall
 * that outerface patch data uses the same indices as face-centered data
 * but the data only exists on the faces that coincide with patch boundaries.
 * It is derived from the hier::TimeInterpolateOperator base class.
 * The interpolation uses FORTRAN numerical routines.
 *
 * The findCoarsenOperator() operator function returns true if the input
 * variable is an outerface double type, and the string is
 * "STD_LINEAR_TIME_INTERPOLATE".
 *
 * @see hier::TimeInterpolateOperator
 */

class OuterfaceDoubleLinearTimeInterpolateOp:
   public hier::TimeInterpolateOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   OuterfaceDoubleLinearTimeInterpolateOp();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~OuterfaceDoubleLinearTimeInterpolateOp();

   /**
    * Return true if the variable and name string match the standard
    * double outerface interpolation; otherwise, return false.
    */
   bool
   findTimeInterpolateOperator(
      const tbox::Pointer<hier::Variable>& var,
      const std::string& op_name) const;

   /**
    * Perform linear time interpolation between two double outerface
    * patch data sources and place result in the destination patch data.
    * Time interpolation is performed on the intersection of the destination
    * patch data and the input box.  The time to which data is interpolated
    * is provided by the destination data.
    */
   void
   timeInterpolate(
      hier::PatchData& dst_data,
      const hier::Box& where,
      const hier::PatchData& src_data_old,
      const hier::PatchData& src_data_new) const;

private:
};

}
}
#endif
