//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-3-0/source/toolbox/database/Serializable.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 2132 $
// Modified:	$LastChangedDate: 2008-04-14 14:51:47 -0700 (Mon, 14 Apr 2008) $
// Description:	An abstract base class for objects than can be serialized
//

#ifndef included_tbox_Serializable
#define included_tbox_Serializable

#include "SAMRAI_config.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"


namespace SAMRAI {
   namespace tbox {


/**
 * Class Serializable is an abstract base class for those objects
 * that can serialize their data to a database.  This class provides 
 * one function that must be implemented in a derived subclass, 
 * putToDatabase(Pointer<Database> db), which should put
 * its data members into the specified database object.
 *
 * Note that the derivation from DescribedClass is virtual.  The
 * reason for this is to avoid dynamic casting problems for smart pointers.  
 * For some objects in SAMRAI, inheritance from Serializable 
 * introduces an additional class hierarchy apart from some other that 
 * may be used to implement the subclass object.  Pointers to base objects
 * may need to be dynamically cast to derived objects in either hierarchy.
 */

class Serializable : public virtual DescribedClass
{
public:
   /**
    * The constructor for the serializable base class does nothing interesting.
    */
   Serializable();

   /**
    * The virtual destructor for the serializable base class does nothing
    * interesting.
    */
   virtual ~Serializable();

   /**
    * This method serializes the object by writing data to the
    * specified database.  
    *
    * NOTE: The asymetry (not having a "getFromDatabase") is from the
    * historical method for doing SAMRAI restart.  The constructor
    * for a Serializable class should get the database to restore
    * state from by making a getRootDatabase call to the
    * RestartManager.
    */
   virtual void putToDatabase(Pointer<Database> db) = 0;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/Serializable.I"
#endif
#endif
