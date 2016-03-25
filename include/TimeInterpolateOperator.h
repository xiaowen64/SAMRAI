//
// File:	TimeInterpolateOperator.h
// Package:	SAMRAI transfer
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Abstract base class for time interpolation operators.
//

#ifndef included_xfer_TimeInterpolateOperator
#define included_xfer_TimeInterpolateOperator

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_hier_Box
#include "Box.h"
#endif
#ifndef included_hier_PatchData
#include "PatchData.h"
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif


namespace SAMRAI {
    namespace xfer {

/**
 * Class TimeInterpolateOperator<DIM> is an abstract base class for each 
 * time interpolation operator used in the SAMRAI framework.  This class
 * defines the interface between numerical interpolation routines and the
 * rest of the framework.  Each concrete time interpolation operator subclass 
 * provide two operations:
 * 


 * - @b (1) an implementation of the time interpolation operation 
 *            appropriate for its corresponding patch data type.
 * - @b (2) a function that determines whether or not the operator 
 *             matches an arbitrary request for a time interpolation operator.
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
 * @see xfer::Geometry
 */

template<int DIM> class TimeInterpolateOperator : public tbox::DescribedClass
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
   virtual ~TimeInterpolateOperator<DIM>();

   /**
    * Return true if the time interpolation operation matches the
    * variable and name string identifier request; false, otherwise.
    */
   virtual bool findTimeInterpolateOperator(
      const tbox::Pointer< hier::Variable<DIM> >& var,
      const string &op_name) const = 0;

   /**
    * Perform time interpolation between two patch data sources 
    * and place result in the destination patch data.  Time interpolation 
    * is performed on the intersection of the destination patch data and 
    * the input box.  The time to which data is interpolated is provided
    * by the destination data.
    */
   virtual void timeInterpolate(
      hier::PatchData<DIM>& dst_data,
      const hier::Box<DIM>& where,
      const hier::PatchData<DIM>& src_data_old,
      const hier::PatchData<DIM>& src_data_new) const = 0;

private:
   // Neither of these is implemented.
   TimeInterpolateOperator(const TimeInterpolateOperator<DIM>&);
   void operator=(const TimeInterpolateOperator<DIM>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "TimeInterpolateOperator.C"
#endif
