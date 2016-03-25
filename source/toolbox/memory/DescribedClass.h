//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/memory/DescribedClass.h $
// Package:	SAMRAI toolbox for RTTI
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1808 $
// Modified:	$LastChangedDate: 2007-12-19 16:38:32 -0800 (Wed, 19 Dec 2007) $
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
 * @brief 
 * Base class for all objects that use run-time type
 * identification (RTTI).
 *
 * This is needed so we can use dynamic casting in our smart pointer
 * implementation.
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
