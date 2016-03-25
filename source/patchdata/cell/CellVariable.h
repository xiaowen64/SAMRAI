//
// File:	CellVariable.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Variable class for defining cell centered variables
//

#ifndef included_pdat_CellVariable
#define included_pdat_CellVariable

#ifndef included_SAMRAI_config
#include "SAMRAI_config.h"
#endif
#ifndef included_tbox_Complex
#include "tbox/Complex.h"
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_hier_Variable
#include "Variable.h"
#endif

namespace SAMRAI {
    namespace pdat {

/**
 * Class CellVariable<DIM> is a templated variable class used to define 
 * cell-centered quantities on an AMR mesh.  It is templated on the type 
 * of the underlying data (e.g., double, int, bool, etc.).  Cell variable 
 * data is associated with the centers of cells.  For example, a 
 * three-dimensional cell variable can be used to create cell-centered
 * data arrays over a box [l0:u0,l1:u1,l2:u2] that can be dimensioned as: 
 * \verbatim

     [ l0 : u0 ,
       l1 : u1 ,
       l2 : u2 , d ]

 * \endverbatim
 * where d is the depth index (i.e., number of values at each cell index 
 * location).  One- and two-dimensional cell variables define storage similarly.
 * For more information on indexing and manipulating cell patch data objects, 
 * see the classes CellData<DIM> and CellGeometry<DIM>.
 *
 * @see pdat::CellData
 * @see pdat::CellDataFactory
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class CellVariable : public hier::Variable<DIM>
{
public:
   /**
    * Create a cell variable object having properties specified by the
    * name and depth (i.e., number of data values at each index location).
    * A default depth of one is provided. 
    */
   CellVariable(const string &name,
                      int depth = 1);

   /**
    * Virtual destructor for cell variable objects.
    */
   virtual ~CellVariable<DIM,TYPE>();

   /**
    * Return true so that the cell data quantities will always be treated as though
    * fine values represent them on coarse-fine interfaces.  Note that this is 
    * really artificial since the cell data index space matches the cell-centered 
    * index space for AMR patches.  Thus, cell data does not live on patch borders
    * and so there is no ambiguity reagrding coarse-fine interface values.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /**
    * Return false since the cell data index space matches the cell-centered
    * index space for AMR patches.  Thus, cell data does not live on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return false;}

private:
   CellVariable(const CellVariable<DIM,TYPE>&);// not implemented
   void operator=(const CellVariable<DIM,TYPE>&);	// not implemented
};


}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "CellVariable.C"
#endif
