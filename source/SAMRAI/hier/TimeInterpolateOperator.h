/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2010 Lawrence Livermore National Security, LLC
 * Description:   Abstract base class for time interpolation operators. 
 *
 ************************************************************************/

#ifndef included_hier_TimeInterpolateOperator
#define included_hier_TimeInterpolateOperator

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/PatchData.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/tbox/DescribedClass.h"
#include "SAMRAI/tbox/Pointer.h"

#include <string>

namespace SAMRAI {
namespace hier {

/**
 * Class TimeInterpolateOperator is an abstract base class for each
 * time interpolation operator used in the SAMRAI framework.  This class
 * defines the interface between numerical interpolation routines and the
 * rest of the framework.  Each concrete time interpolation operator subclass
 * provide two operations:
 *
 *
 *
 * - @b (1) an implementation of the time interpolation operation
 *            appropriate for its corresponding patch data type.
 * - @b (2) a function that determines whether or not the operator
 *             matches an arbitrary request for a time interpolation operator.
 *
 *
 *
 * To add a new time interpolation operator (either for a new patch data type
 * or for a new time interpolation routine on an existing type), define
 * the operator by inheriting from this abstract base class.  The operator
 * subclass must implement the interpolation operation in the timeInterpolate()
 * function, and provide a response to a general operator request in the
 * findTimeInterpolateOperator() function.  Then, the new operator must be
 * added to the operator list for the appropriate transfer geometry object
 * using the Geometry<DIM>::addTimeInterpolateOperator() function.
 *
 * Although time interpolation operators usually depend only on patch data
 * centering and data type and not the mesh coordinate system, they are
 * defined in the @em geometry package.
 *
 * @see hier::TransferOperatorRegistry
 */

class TimeInterpolateOperator:public tbox::DescribedClass
{
public:
   /**
    * The default constructor for the coarsening operator does
    * nothing interesting.
    */
   TimeInterpolateOperator();

   /**
    * The virtual destructor for the coarsening operator does
    * nothing interesting.
    */
   virtual ~TimeInterpolateOperator();

   /**
    * Return true if the time interpolation operation matches the
    * variable and name std::string identifier request; false, otherwise.
    */
   virtual bool
   findTimeInterpolateOperator(
      const tbox::Pointer<Variable>& var,
      const std::string& op_name) const = 0;

   /**
    * Perform time interpolation between two patch data sources
    * and place result in the destination patch data.  Time interpolation
    * is performed on the intersection of the destination patch data and
    * the input box.  The time to which data is interpolated is provided
    * by the destination data.
    */
   virtual void
   timeInterpolate(
      PatchData& dst_data,
      const Box& where,
      const PatchData& src_data_old,
      const PatchData& src_data_new) const = 0;

private:
   // Neither of these is implemented.
   TimeInterpolateOperator(
      const TimeInterpolateOperator&);
   void
   operator = (
      const TimeInterpolateOperator&);

};

}
}
#endif
