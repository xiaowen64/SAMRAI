//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/mathops/array/ArrayDataNormOpsInteger.h $
// Package:	SAMRAI mathops
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Templated norm operations for real data arrays.
//

#ifndef included_math_ArrayDataNormOpsInteger
#define included_math_ArrayDataNormOpsInteger

#include "SAMRAI_config.h"
#include "Box.h"
#include "ArrayData.h"

namespace SAMRAI {
    namespace math {

/**
 * Class ArrayDataNormOpsInteger<DIM> provides a set of common norm 
 * operations that may be applied to arrays of integer data values
 * maintained as pdat::ArrayData<DIM> objects.  The intent of this class 
 * is to provide a single implementation of these operations as they are needed 
 * by objects that perform these operations on the standard array-based patch 
 * data types (i.e., cell-centered, face-centered, node-centered).  
 * Note that each operation is performed on the intersection of the box in 
 * the function argument list and the boxes associated with all 
 * pdat::ArrayData<DIM> objects.  Currently, the only norm operation implemented 
 * in this class is the absolute value operation.
 *
 * @see pdat::ArrayData
 */

template<int DIM>
class ArrayDataNormOpsInteger
{
public:
   /** 
    * Empty constructor and destructor.
    */
   ArrayDataNormOpsInteger();

   ~ArrayDataNormOpsInteger<DIM>();

   /**
    * Set destination component to absolute value of source component.  
    * That is, each destination entry is set to \f$d_i = \| s_i \|\f$.
    */
   void abs(pdat::ArrayData<DIM,int>& dst,
            const pdat::ArrayData<DIM,int>& src,
            const hier::Box<DIM>& box) const;

private:
   // The following are not implemented:
   ArrayDataNormOpsInteger(const ArrayDataNormOpsInteger<DIM>&);
   void operator=(const ArrayDataNormOpsInteger<DIM>&);
};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "ArrayDataNormOpsInteger.C"
#endif
