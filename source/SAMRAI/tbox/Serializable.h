/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and COPYING.LESSER.
 *
 * Copyright:     (c) 1997-2011 Lawrence Livermore National Security, LLC
 * Description:   An abstract base class for objects to be serialized
 *
 ************************************************************************/

#ifndef included_tbox_Serializable
#define included_tbox_Serializable

#include "SAMRAI/SAMRAI_config.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/DescribedClass.h"

#include <boost/shared_ptr.hpp>

namespace SAMRAI {
namespace tbox {

/*!
 * @brief Class Serializable is an abstract base class to provide
 * database serialization for objects.
 *
 * Objects needing to serialize their data to a database implement the
 * function putToDatabase(boost::shared_ptr<tbox::Database>) specifying
 * how to put its data members into the database instance.
 *
 * @note The interface asymmetric in that is no corresponding
 * "getFromDatabase" function.  This stems from SAMRAI's historical
 * restart procedure and object creation is object initialization
 * pattern.  Thus restarting is done in the constructors for classes
 * implementing Serializable not by an additional method.  The
 * database used to restore state from is obtained via a
 * getRootDatabase() call to the RestartManager.
 *
 * @note The derivation from DescribedClass is virtual.  The reason
 * for this is to avoid dynamic casting problems for smart pointers.
 * For some objects in SAMRAI, inheritance from Serializable
 * introduces a class hierarchy that is distinct from the class
 * hierarchy used to implement the class.  Pointers to base objects
 * need to be dynamically cast to derived objects in either hierarchy.
 */

class Serializable:public virtual DescribedClass
{
public:
   /*!
    * Default constructor.
    */
   SAMRAI_INLINE_KEYWORD Serializable() {
   }

   /*!
    * Virtual destructor
    */
   virtual ~Serializable() {
   }

   /*!
    * @brief Serializes the object by writing necessary state to the
    * specified database.
    *
    */
   virtual void
   putToDatabase(
      boost::shared_ptr<Database> database) = 0;

};

}
}

#endif
