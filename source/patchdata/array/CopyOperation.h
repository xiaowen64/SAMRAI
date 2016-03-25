//
// File:	$URL$
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	Copy operation on single array data elements templated on data type
//

#ifndef included_pdat_CopyOperation
#define included_pdat_CopyOperation

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

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

#endif
