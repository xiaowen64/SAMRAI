/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright 
 * information, see COPYRIGHT and COPYING.LESSER. 
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   Base class for run-time type identification 
 *
 ************************************************************************/

#ifndef included_tbox_DescribedClass
#define included_tbox_DescribedClass

#include "SAMRAI/SAMRAI_config.h"

namespace SAMRAI {
namespace tbox {

/**
 * @brief
 * Required base class for all classes to be used with smart pointers.
 *
 * All of the classes pointed to by a SAMRAI smart pointer need to be
 * derived from this class.  SAMRAI smart pointers require that the
 * C++ run-time type identification (RTTI) dynamic casting operator (
 * dynamic_cast() ) work on the objects being pointed to.  The dynamic
 * cast operator requires the classes are in the same class hierarchy.
 * Deriving from DesribedClasss guarantees this hence the reason for
 * this unattractive requirement.
 *
 * Notes:
 *
 * The SAMRAI developers recognize this requirement is a not ideal and
 * that the dynamic_cast that is being done implicitly in
 * the smart pointer assignment is considered to be a poor choice by
 * many other reference counting pointer implementations.  With
 * hindsight a different approach would have been used, likely more
 * consistent with the Boost reference counting pointers.
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
