//
// File:	DescribedClass.h
// Package:	SAMRAI toolbox for RTTI
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	Base class for run-time type identification
//

#ifndef included_tbox_DescribedClass
#define included_tbox_DescribedClass

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif

namespace SAMRAI {
   namespace tbox {

/**
 * Class DescribedClass is a base class for all objects that use
 * run-time type identification (RTTI).  This is needed so we can use
 * dynamic casting in our smart pointer implementation.
 *
 * @see tbox::Pointer
 */

class DescribedClass
{
public:
   /**
    * The virtual destructor for DescribedClass does nothing interesting.
    */
   virtual ~DescribedClass();
};


}
}

#endif
