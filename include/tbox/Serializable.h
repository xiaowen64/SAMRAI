//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-0/source/toolbox/database/Serializable.h $
// Package:	SAMRAI toolbox
// Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
// Revision:	$LastChangedRevision: 1704 $
// Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
// Description:	An abstract base class for objects than can be serialized
//

#ifndef included_tbox_Serializable
#define included_tbox_Serializable

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_DescribedClass
#include "tbox/DescribedClass.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif


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
    */
   virtual void putToDatabase(Pointer<Database> db) = 0;

};

}
}

#ifndef DEBUG_NO_INLINE
#include "tbox/Serializable.I"
#endif
#endif
