//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/memory/DescribedClass.h $
// Package:	SAMRAI toolbox for RTTI
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	Base class for run-time type identification
//

#ifndef included_tbox_DescribedClass
#define included_tbox_DescribedClass

#include "SAMRAI_config.h"

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
