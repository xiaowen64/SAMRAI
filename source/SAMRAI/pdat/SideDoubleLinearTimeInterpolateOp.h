/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Linear time interp operator for side-centered double patch data. 
 *
 ************************************************************************/

#ifndef included_pdat_SideDoubleLinearTimeInterpolateOp
#define included_pdat_SideDoubleLinearTimeInterpolateOp

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/TimeInterpolateOperator.h"
#include "SAMRAI/tbox/Pointer.h"

#include <string>

namespace SAMRAI {
namespace pdat {

/**
 * Class SideDoubleLinearTimeInterpolateOp implements standard
 * linear time interpolation for side-centered double patch data.
 * It is derived from the hier::TimeInterpolateOperator base class.
 * The interpolation uses FORTRAN numerical routines.
 *
 * The findCoarsenOperator() operator function returns true if the input
 * variable is a side-centered double, and the string is
 * "STD_LINEAR_TIME_INTERPOLATE".
 *
 * @see hier::TimeInterpolateOperator
 */

class SideDoubleLinearTimeInterpolateOp:
   public hier::TimeInterpolateOperator
{
public:
   /**
    * Uninteresting default constructor.
    */
   SideDoubleLinearTimeInterpolateOp();

   /**
    * Uninteresting virtual destructor.
    */
   virtual ~SideDoubleLinearTimeInterpolateOp();

   /**
    * Return true if the variable and name string match the standard
    * side-centered double interpolation; otherwise, return false.
    */
   bool
   findTimeInterpolateOperator(
      const tbox::Pointer<hier::Variable>& var,
      const std::string& op_name) const;

   /**
    * Perform linear time interpolation between two side-centered double
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
