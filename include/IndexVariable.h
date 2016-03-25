//
// File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-4-0/source/patchdata/index/IndexVariable.h $
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2008 Lawrence Livermore National Security, LLC
// Release:	0.1
// Revision:	$LastChangedRevision: 2224 $
// Modified:	$LastChangedDate: 2008-06-20 17:51:16 -0700 (Fri, 20 Jun 2008) $
// Description:	hier::Variable class for defining irregular index variables
//

#ifndef included_pdat_IndexVariable
#define included_pdat_IndexVariable

#include "SAMRAI_config.h"
#include "IntVector.h"
#include "Variable.h"
#ifndef included_String
#include <string>
#define included_String
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class IndexVariable<DIM,TYPE,BOX_GEOMETRY> is a templated variable
 * class used to define quantities that exist on an irregular
 * cell-centered index set.  The template parameter TYPE defines the
 * storage at each index location.  For example, this class is used to
 * represent embedded boundary features as a regular patch data type
 * using the BoundaryCell class as the template type.  The template
 * parameter BOX_GEOMETRY allows IndexVariables to be instantiated
 * with a provided centering and geometry in index space via a
 * BoxGeometry (e.g. CellGeometry, NodeGeometry).
 *
 * Please consult the README file in the index data source directory for 
 * instructions on using this class to provide other irregular index set
 * types.
 *
 * @see pdat::IndexData
 * @see pdat::IndexDataFactory
 * @see pdat::Variable
 */

template<int DIM, class TYPE, class BOX_GEOMETRY>
class IndexVariable : public hier::Variable<DIM>
{
public:
   /**
    * Create an index variable object with the specified name.
    */
   IndexVariable(const std::string &name);

   /**
    * Virtual destructor for index variable objects.
    */
   virtual ~IndexVariable<DIM,TYPE,BOX_GEOMETRY>();

   /**
    * Return true so that the index data quantities will always be treated as cell-
    * centered quantities as far as communication is concerned.  Note that this is
    * really artificial since the cell data index space matches the cell-centered
    * index space for AMR patches.  Thus, cell data does not live on patch borders
    * and so there is no ambiguity reagrding coarse-fine interface values.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /**
    * Return false since the index data index space matches the cell-centered
    * index space for AMR patches.  Thus, index data does not live on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return false;}

private:
   IndexVariable(const IndexVariable<DIM,TYPE,BOX_GEOMETRY>&); // not implemented
   void operator=(const IndexVariable<DIM,TYPE,BOX_GEOMETRY>&);      // not implemented

};



}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "IndexVariable.C"
#endif
