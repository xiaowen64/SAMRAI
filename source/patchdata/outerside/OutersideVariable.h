//
// File:	OutersideVariable.h
// Package:	SAMRAI patch data
// Copyright:	(c) 1997-2005 The Regents of the University of California
// Revision:	$Revision: 173 $
// Modified:	$Date: 2005-01-19 09:09:04 -0800 (Wed, 19 Jan 2005) $
// Description:	hier::Variable class for defining outerside centered variables
//

#ifndef included_pdat_OutersideVariable
#define included_pdat_OutersideVariable

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
 * Class OutersideVariable<DIM> is a templated variable class used to define
 * side-centered quantities on an AMR mesh.  It is templated on the type
 * of the underlying data (e.g., double, int, bool, etc.).  Outerside variable
 * data is associated with the sides (or faces) of cells.  However, it differs
 * from the SideVariable<DIM> class in that outerside quantities reside only
 * on the sides residing on the boundary of a patch.  Outerside data is
 * stored in 2*DIM arrays, each of which holds values for sides having the
 * same outward normal vector.  For example, a three-dimensional outerside
 * variable can be used to create side-centered data arrays over a box
 * [l0:u0,l1:u1,l2:u2] that can be dimensioned as:
 * \verbatim

     [ l1 : u1 ,
       l2 : u2 , d ]   ,

     [ l0 : u0 ,
       l2 : u2 , d ]   ,

     [ l0 : u0 ,
       l1 : u1 , d ]   ,

 * \endverbatim
 * for the upper and lower x, y, and z (or 0, 1, 2) face directions,
 * respectively, and where d is the depth index (i.e., number of values at
 * each side index location).  Note that the array orderings match the 
 * conventions of the SideVariable<DIM> class.  One- and two-dimensional 
 * outerside variables define storage similarly.  For more information on 
 * indexing and manipulating outerside patch data objects, see the classes 
 * OutersideData<DIM> and OutersideGeometry<DIM>.
 *
 * IMPORTANT: The class OuterfaceVariable<DIM> and associated classes define
 * the same storage as this outerside variable class, except that the indices
 * are not permuted in the outerside data type.  Also, outerface and outerside
 * data classes are intended to interact with their face-centered and
 * side-centered data counterparts, respectively.  Mixing types, while
 * allowed, is discouraged to prevent undesirable behavior.
 *
 * @see pdat::SideData
 * @see pdat::OutersideData
 * @see pdat::OutersideDataFactory
 * @see hier::Variable
 */

template<int DIM, class TYPE>
class OutersideVariable : public hier::Variable<DIM>
{
public:
   /**
    * Create an outerside variable object having properties specified by the
    * name and depth (i.e., number of data values at each index location).
    * The default depth is one.   The ghost cell width for all outerside
    * data is currently fixed at zero; this may be changed in the future
    * if needed.
    */
   OutersideVariable(const string &name, 
                           int depth = 1);

   /**
    * Virtual destructor for outerside variable objects.
    */
   virtual ~OutersideVariable<DIM,TYPE>();

   /**
    * Return a boolean true value indicating that fine patch values take precedence
    * on coarse-fine interfaces.
    */
   bool fineBoundaryRepresentsVariable() const {return true;}

   /**
    * Return true since the side data index space (and hence the outerside data index
    * space) extends beyond the interior of patches.  That is, outerside data lives
    * on patch borders.
    */
   bool dataLivesOnPatchBorder() const {return true;}

private:
   // neither of the following functions are implemented
   OutersideVariable(const OutersideVariable<DIM,TYPE>&);
   void operator=(const OutersideVariable<DIM,TYPE>&);

};

}
}
#endif

#ifdef INCLUDE_TEMPLATE_IMPLEMENTATION
#include "OutersideVariable.C"
#endif
