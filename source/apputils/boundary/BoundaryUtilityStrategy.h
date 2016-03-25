//
// File:        BoundaryUtilityStrategy
// Package:     SAMRAI application utilities
// Copyright:   (c) 1997-2005 The Regents of the University of California
// Revision:    $Revision: 173 $
// Modified:    $Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description: Interface for processing user-defined boundary data in
//              CartesianBoundaryUtilities classes
//

#ifndef included_appu_BoundaryUtilityStrategy
#define included_appu_BoundaryUtilityStrategy

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Database
#include "tbox/Database.h"
#endif
#ifndef included_tbox_Pointer
#include "tbox/Pointer.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_tbox_Utilities
#include "tbox/Utilities.h"
#endif

namespace SAMRAI {
    namespace appu {

/*!
 * @brief Class BoundaryUtilityStrategy is an abstract base 
 * class that declares an interface that allows application code to
 * read problem-specific boundary data when using the SAMRAI boundary 
 * utilities.  Currently, there are two virtual member functions defined.  
 * One allows users to read problem-specific DIRICHLET boundary values 
 * from an input database; the other does the same for NEUMANN boundary 
 * values.  More virtual functions may be added in the future 
 * as additional boundary conditions are supported.
 *
 * See the include file BoundaryDefines.h for integer constant
 * definitions that apply for the various boundary types, locations,
 * and boundary conditions.
 *
 * @see appu::CartesianBoundaryUtilities2
 * @see appu::CartesianBoundaryUtilities3
 */

class BoundaryUtilityStrategy 
{
public:
   /*!
    * The default constructor for BoundaryUtilityStrategy does nothing
    * interesting.
    */
   BoundaryUtilityStrategy();

   /*!
    * The destructor for BoundaryUtilityStrategy does nothing
    * interesting.
    */
   virtual ~BoundaryUtilityStrategy();

   /*!
    * Read DIRICHLET boundary state values for an edge (in 2d) or a face
    * (in 3d) from a given database.  This virtual function is given a 
    * blank implementation here to avoid the need for users to do the same 
    * if they do not need this functionality. 
    * 
    * @param db      Input database from which to read boundary values.
    * @param db_name Name of input database (e.g., for error reporting).
    * @param bdry_location_index Integer index for location of edge (in 2d)
    *                            or face (in 3d) boundary.
    */
   virtual void readDirichletBoundaryDataEntry(
      tbox::Pointer<tbox::Database> db,
      string& db_name,
      int bdry_location_index)
   {
      NULL_USE(db);
      NULL_USE(db_name);
      NULL_USE(bdry_location_index);
   } 

   /*!
    * Read NEUMANN boundary state values for an edge (in 2d) or a face
    * (in 3d) from a given database.  This virtual function is given a
    * blank implementation here to avoid the need for users to do the same
    * if they do not need this functionality.
    *
    * @param db      Input database from which to read boundary values.
    * @param db_name Name of input database (e.g., for error reporting).
    * @param bdry_location_index Integer index for location of edge (in 2d)
    *                            or face (in 3d) boundary.
    */
   virtual void readNeumannBoundaryDataEntry(
      tbox::Pointer<tbox::Database> db,
      string& db_name,
      int bdry_location_index)
   {
      NULL_USE(db);
      NULL_USE(db_name);
      NULL_USE(bdry_location_index);
   }

};

}
}
#endif

