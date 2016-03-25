//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	Sum operation on single array data elements templated on data type
//

#ifndef included_pdat_SumOperation
#define included_pdat_SumOperation

#include "SAMRAI_config.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class SumOperation<TYPE> encapsulates a summation 
 * operation into an object.
 */

template <class TYPE>
class SumOperation
{
public:
   /*!
    * The default constructor does nothing interesting.
    */
   SumOperation();

   /*!
    * The destructor does nothing interesting.
    */
   ~SumOperation();

   /*!
    * The operator adds the source value to the destination.
    */
   void operator()(TYPE& vdst, const TYPE& vsrc) const;

private:
   SumOperation(const SumOperation<TYPE>&);   // not implemented
   void operator=(const SumOperation<TYPE>&);  // not implemented
};

}
}

#ifndef DEBUG_NO_INLINE
#include "SumOperation.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "SumOperation.C"
#endif

#endif
