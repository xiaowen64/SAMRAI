//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Sum operation on single array data elements templated on data type
//

#ifndef included_pdat_SumOperation
#define included_pdat_SumOperation

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

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

#endif
