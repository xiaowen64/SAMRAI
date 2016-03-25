//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2195 $
// Modified:	$LastChangedDate: 2008-05-14 11:33:30 -0700 (Wed, 14 May 2008) $
// Description:	Copy operation on single array data elements templated on data type
//

#ifndef included_pdat_CopyOperation
#define included_pdat_CopyOperation

#include "SAMRAI_config.h"

namespace SAMRAI {
    namespace pdat {

/*!
 * @brief Class CopyOperation<TYPE> encapsulates a copy
 * operation into an object.
 */

template <class TYPE>
class CopyOperation
{
public:
   /*!
    * The default constructor does nothing interesting.
    */
   CopyOperation();

   /*!
    * The destructor does nothing interesting.
    */
   ~CopyOperation();

   /*!
    * The operator copies the source value to the destination.
    */
   void operator()(TYPE& vdst, const TYPE& vsrc) const;

private:
   CopyOperation(const CopyOperation<TYPE>&);	// not implemented
   void operator=(const CopyOperation<TYPE>&);	// not implemented
};


}
}

#ifndef DEBUG_NO_INLINE
#include "CopyOperation.I"
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CopyOperation.C"
#endif

#endif
